#include "rand_flash.hpp"
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
using namespace randflash;
using namespace ls;

// ----------------------------------------------------------------------------------
// 完善后的三相含水体系初始化：顺序拆分策略 (Sequential Splitting)
// 1. 识别水组分。
// 2. 假设水相主要由水组成，烃相由两相 RR 估算气液比。
// 3. 构建合成 K 值，通过严格的分配公式生成满足元素守恒的初始 n。
// ----------------------------------------------------------------------------------
InitResult RandFlash::initializeThreePhaseWater(
    const SystemContext& sys,
    double margin) const
{
    const size_t C = sys.feedMoles.size();
    // 计算总进料量
    const double Ftot = std::accumulate(sys.feedMoles.begin(), sys.feedMoles.end(), 0.0);
    
    // 计算归一化进料 z
    std::vector<double> z = sys.feedMoles;
    if (Ftot > 1e-20) {
        for(double& val : z) val /= Ftot;
    }

    // ====================================================
    // 1. 识别水组分索引
    // ====================================================
    auto names = thermo_.getComponentNames();
    int water_idx = -1;
    for(size_t i=0; i<C; ++i) {
        std::string n = names[i];
        std::transform(n.begin(), n.end(), n.begin(), ::toupper);
        if(n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
            water_idx = (int)i;
            break;
        }
    }

    // ====================================================
    // 2. 估算 K 值 (基准相：油相 Oil)
    //    K_gas = y / x_oil (使用 Wilson 公式)
    //    K_aq  = x_aq / x_oil (构造不混溶 K 值)
    // ====================================================
    std::vector<double> K_gas(C), K_aq(C);
    
    // 2.1 计算 Wilson K 值 (用于气-油平衡)
    thermo_.wilsonK(sys.temperature, sys.pressure, K_gas);

    // 2.2 构造水-油 K 值 (极不混溶假设)
    // 经验值：水在水相/油相的分配系数极大，烃类极小
    const double K_w_in_aq = 1.0e4;  // 水组分倾向于水相
    const double K_hc_in_aq = 1.0e-5; // 烃类组分极难溶于水

    for(size_t i=0; i<C; ++i) {
        if ((int)i == water_idx) {
            K_aq[i] = K_w_in_aq; 
            // 修正气相 K 值：水在气相也比在油相多(蒸汽压)，但不如水相多
            // 这里保留 Wilson 计算的值即可，通常 Wilson 对水的预测在低压下尚可
        } else {
            K_aq[i] = K_hc_in_aq;
        }
    }

    // ====================================================
    // 3. 估算相分率 Beta (顺序拆分法)
    // ====================================================
    double betaW = 0.0;
    double betaV = 0.0;
    double betaO = 0.0;

    // 3.1 估算水相分率 (假设所有水组分都进入水相)
    if (water_idx >= 0) {
        betaW = z[water_idx]; 
    }
    
    // 边界保护：如果全是水或没水
    if (betaW > 0.999) betaW = 0.999;
    if (betaW < 1e-6)  betaW = 1e-6; // 即使没水，也给一点点，防止矩阵奇异

    // 3.2 对剩余烃类部分进行 2相 Rachford-Rice 估算气/油比
    // 剩余烃类总分率
    double z_hc_total = 1.0 - betaW;
    
    // 归一化的烃类组成 z_hc
    std::vector<double> z_hc(C, 0.0);
    bool has_hc = false;
    for(size_t i=0; i<C; ++i) {
        if ((int)i != water_idx) {
            z_hc[i] = z[i] / z_hc_total;
            if (z[i] > 1e-10) has_hc = true;
        }
    }

    double alpha_vap = 0.5; // 烃类部分的气相分率 (V / (V+O))

    if (has_hc) {
        // 求解简单的 Rachford-Rice: sum (zi * (K-1) / (1 + alpha*(K-1))) = 0
        // 仅针对非水组分
        auto rr_func = [&](double a) {
            double sum = 0.0;
            for(size_t i=0; i<C; ++i) {
                if ((int)i == water_idx) continue;
                double denom = 1.0 + a * (K_gas[i] - 1.0);
                sum += z_hc[i] * (K_gas[i] - 1.0) / denom;
            }
            return sum;
        };

        // 二分法求解 alpha_vap
        double low = 0.0, high = 1.0;
        double f_low = rr_func(low);
        double f_high = rr_func(high);

        if (f_low * f_high < 0.0) {
            for(int k=0; k<20; ++k) {
                double mid = 0.5 * (low + high);
                double f_mid = rr_func(mid);
                if (std::abs(f_mid) < 1e-6) { alpha_vap = mid; break; }
                if (f_mid * f_low > 0.0) { low = mid; f_low = f_mid; }
                else { high = mid; }
            }
            alpha_vap = 0.5 * (low + high);
        } else {
            // 单相区：根据函数值符号判断是过冷液还是过热气
            // f_low > 0 说明全是气 (bubble point < P) -> alpha = 1.0 (近似)
            // f_high < 0 说明全是液 -> alpha = 0.0 (近似)
            if (f_low > 0) alpha_vap = 0.99;
            else alpha_vap = 0.01;
        }
    }

    // 3.3 转换为全局相分率
    betaV = alpha_vap * (1.0 - betaW);
    betaO = 1.0 - betaV - betaW;

    // 最小相分率保护
    if (betaV < 1e-4) betaV = 1e-4;
    if (betaO < 1e-4) betaO = 1e-4;
    
    // 重新归一化 Beta
    double sum_beta = betaV + betaO + betaW;
    betaV /= sum_beta;
    betaO /= sum_beta;
    betaW /= sum_beta;

    // ====================================================
    // 4. 生成 n (严格满足元素守恒)
    //    n_{i,j} = Ftot * z_i * (beta_j * K_{ij}) / (sum_k beta_k K_{ik})
    // ====================================================
    InitResult out;
    out.n_phases.resize(3);
    out.n_phases[0].resize(C); // Vapor (相索引 0)
    out.n_phases[1].resize(C); // Oil   (相索引 1)
    out.n_phases[2].resize(C); // Aqueous (相索引 2)
    out.beta.resize(3);
    out.compositions.resize(3);

    // K 值定义回顾：
    // K_gas[i] = y[i] / x_oil[i]  =>  K_{i,0} / K_{i,1}
    // K_aq[i]  = x_aq[i] / x_oil[i] => K_{i,2} / K_{i,1}
    // 令 K_{i,Oil} = 1.0, 则 K_{i,Vapor} = K_gas[i], K_{i,Aq} = K_aq[i]

    for(size_t i=0; i<C; ++i) {
        // 分母 D_i = sum_j (beta_j * K_{ij})
        double Di = betaO * 1.0 + betaV * K_gas[i] + betaW * K_aq[i];
        
        // 防止除零
        if (Di < 1e-30) Di = 1e-30;

        double total_moles_i = sys.feedMoles[i]; // z_i * Ftot

        // Oil 相摩尔数 (基准)
        double n_oil = total_moles_i * (betaO * 1.0) / Di;

        out.n_phases[1][i] = n_oil;                     // Oil
        out.n_phases[0][i] = n_oil * (betaV * K_gas[i] / betaO); // Vapor
        out.n_phases[2][i] = n_oil * (betaW * K_aq[i] / betaO);  // Aqueous
    }

    // ====================================================
    // 5. 强制全局缩放 (Final Re-scaling)
    //    消除浮点数累加误差，确保 residual < 1e-15
    // ====================================================
    for(size_t i=0; i<C; ++i) {
        double current_sum = out.n_phases[0][i] + out.n_phases[1][i] + out.n_phases[2][i];
        if (current_sum > 1e-20) {
            double scale = sys.feedMoles[i] / current_sum;
            out.n_phases[0][i] *= scale;
            out.n_phases[1][i] *= scale;
            out.n_phases[2][i] *= scale;
        }
    }

    // 6. 填充结果元数据 (beta, x)
    for(int j=0; j<3; ++j) {
        out.beta[j] = std::accumulate(out.n_phases[j].begin(), out.n_phases[j].end(), 0.0);
        out.compositions[j] = out.n_phases[j];
        if (out.beta[j] > 1e-20) {
            double invBeta = 1.0 / out.beta[j];
            for(double& val : out.compositions[j]) val *= invBeta;
        } else {
            std::fill(out.compositions[j].begin(), out.compositions[j].end(), 1.0/C);
        }
    }

    // 打印调试信息
    std::cout << "[Init] Sequential Split Init:\n"
              << "       Est. Beta V: " << betaV << " (Alpha_vap of HC: " << alpha_vap << ")\n"
              << "       Est. Beta O: " << betaO << "\n"
              << "       Est. Beta W: " << betaW << " (Water idx: " << water_idx << ")\n";

    return out;
}

  MultiFlashResult RandFlash::solveMultiPhaseCore(
  const FlashInput& input,
  int nPhases,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<std::vector<double>>& initialPhaseCompositions,
  int maxIter,
  double tol)
{
    // 1. 构建 Context
    SystemContext sys;
    sys.temperature = input.temperature;
    sys.pressure = input.pressure;
    sys.feedMoles = input.feedMoles;
    sys.elementMatrix = elementMatrix;
    sys.phases.resize(nPhases);

    const size_t C = input.feedMoles.size();
    bool initialized = false;

    // 情况 A: 两相且无初值 -> 调用标准两相初始化
    if (nPhases == 2 && initialPhaseCompositions.empty()) {
        auto init = initializeTwoPhase(sys, {}, {});
        sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], thermo_.vaporPhaseFlag()};
        sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], thermo_.liquidPhaseFlag()};
        initialized = true;
    } 
    // 情况 B: 三相且无初值 -> 调用含水三相特殊初始化 (需配合之前给出的 initializeThreePhaseWater 实现)
    else if (nPhases == 3 && initialPhaseCompositions.empty()) {
        // 确保你的类中已经添加了 initializeThreePhaseWater
        auto init = initializeThreePhaseWater(sys);
        sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], thermo_.vaporPhaseFlag()};
        sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], thermo_.liquidPhaseFlag()};
        sys.phases[2].state = {input.temperature, input.pressure, init.n_phases[2], thermo_.liquidPhaseFlag()};
        initialized = true;
    }
    // 情况 C: 有用户提供的初值 -> 最小二乘投影 + 强制守恒缩放
    else if (!initialPhaseCompositions.empty()) {
        if ((int)initialPhaseCompositions.size() != nPhases) {
             std::cerr << "[Warning] Initial compositions count != nPhases. Using default logic.\n";
        } else {
            // 1. 准备归一化的组成矩阵 X (C x F)
            std::vector<std::vector<double>> X_norm(nPhases, std::vector<double>(C));
            for(int j=0; j<nPhases; ++j) {
                double s = std::accumulate(initialPhaseCompositions[j].begin(), initialPhaseCompositions[j].end(), 0.0);
                if (s < 1e-15) s = 1.0; 
                for(size_t i=0; i<C; ++i) X_norm[j][i] = initialPhaseCompositions[j][i] / s;
            }

            // 2. 构建正规方程 A * beta = b
            // A = X^T * X (F x F), b = X^T * z (F x 1)
            // 这里的 z 是 input.feedMoles
            std::vector<double> A_ls(nPhases * nPhases, 0.0);
            std::vector<double> b_ls(nPhases, 0.0);

            for (int j = 0; j < nPhases; ++j) {
                // 构建 b_ls[j] = x_j . z
                for (size_t i = 0; i < C; ++i) {
                    b_ls[j] += X_norm[j][i] * input.feedMoles[i];
                }
                // 构建 A_ls[j][k] = x_j . x_k
                for (int k = 0; k < nPhases; ++k) {
                    double dot = 0.0;
                    for (size_t i = 0; i < C; ++i) {
                        dot += X_norm[j][i] * X_norm[k][i];
                    }
                    A_ls[j * nPhases + k] = dot;
                }
            }

            // 3. 求解 Beta
            // 使用传入的线性求解器解小规模方程
            double dummy_res = 0.0;
            std::vector<double> betas = linearSolver_.solveDense(nPhases, A_ls, b_ls, &dummy_res);

            // 4. Beta 修正 (防止负值或过小)
            double total_feed = std::accumulate(input.feedMoles.begin(), input.feedMoles.end(), 0.0);
            for(double& b : betas) {
                if(b < 1e-5 * total_feed) b = 1e-5 * total_feed;
            }

            // 5. 生成初步 n 并进行强制行缩放 (Row-wise Scaling)
            // 这是保证 Element Balance 残差为 0 的关键步骤
            std::vector<std::vector<double>> n_temp(nPhases, std::vector<double>(C));
            std::vector<double> n_sum(C, 0.0);

            // 5.1 初步计算 n = beta * x
            for(int j=0; j<nPhases; ++j) {
                for(size_t i=0; i<C; ++i) {
                    n_temp[j][i] = betas[j] * X_norm[j][i];
                    n_sum[i] += n_temp[j][i];
                }
            }

            // 5.2 强制缩放：n_final = n_temp * (z_i / sum_n_i)
            for(size_t i=0; i<C; ++i) {
                double target = input.feedMoles[i];
                double current = n_sum[i];
                double scale = (current > 1e-20) ? (target / current) : 0.0;
                
                // 如果当前组分在所有相中都为0，但进料不为0，则平均分配
                if (current <= 1e-20 && target > 1e-20) {
                    for(int j=0; j<nPhases; ++j) n_temp[j][i] = target / nPhases;
                } else {
                    for(int j=0; j<nPhases; ++j) n_temp[j][i] *= scale;
                }
            }

            // 6. 赋值给 SystemContext
            for(int j=0; j<nPhases; ++j) {
                int flag = (j==0 ? thermo_.vaporPhaseFlag() : thermo_.liquidPhaseFlag());
                sys.phases[j].state = {input.temperature, input.pressure, n_temp[j], flag};
            }
            initialized = true;
            std::cout << "[Init] Initialized " << nPhases << " phases with Least-Squares + Strict Element Scaling.\n";
        }
    }

    // 情况 D: 兜底逻辑
    if (!initialized) {
        std::cout << "[Init] Fallback: Evenly splitting feed.\n";
        std::vector<double> n = input.feedMoles;
        double split = 1.0 / nPhases;
        for(auto& val : n) val *= split;
        for(int j=0; j<nPhases; ++j) {
            int flag = (j==0 ? thermo_.vaporPhaseFlag() : thermo_.liquidPhaseFlag());
            sys.phases[j].state = {input.temperature, input.pressure, n, flag};
        }
    }

    // 3. 核心求解
    return solveGeneral(sys, maxIter, tol);
}

MultiFlashResult RandFlash::solveMultiPhase(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions,
    int maxIter,
    double tol)
{
  // 目前：如果用户没给初值，就默认 2 相；
  // 如果给了，就用给的相数；但核心仍只接受 2 相。
  int F = initialPhaseCompositions.empty()
        ? 2
        : static_cast<int>(initialPhaseCompositions.size());

  return solveMultiPhaseCore(
      input,
      F,
      elementMatrix,
      initialPhaseCompositions,
      maxIter,
      tol);
}

MultiFlashResult RandFlash::solveGeneral(
    SystemContext& sys,
    int maxIter,
    double tol)
{
    MultiFlashResult result;
    result.pressure = sys.pressure;
    result.temperature = sys.temperature;
    result.feedComposition = sys.feedMoles;
    
    const size_t F = sys.phases.size();
    const size_t E = sys.elementMatrix.size();

    std::vector<double> Acoef, rhs, sol;
    std::vector<double> Lambda;
    std::vector<double> deltaBeta(F);
    std::vector<std::vector<double>> dnPhases(F);

    std::cout << "进入通用迭代循环 (F=" << F << "), maxIter = " << maxIter << std::endl;

    try {
        for (int iter = 0; iter < maxIter; ++iter) {
            std::cout << "\n=== RandFlash General Iter " << (iter + 1) << " ===\n";

            std::vector<std::vector<std::vector<double>>> Ms(F);
            std::vector<std::vector<double>> mus(F);
            std::vector<std::vector<double>> nPhases(F);

            for(size_t j=0; j<F; ++j) {
                updatePhaseChemistry(sys.phases[j]);
                
                // 【核心修复】提高特征值下限阈值
                // 1e-10 -> 1e-3
                // 这限制了不稳定相的修正步长幅度，防止产生过大的 delta_n 导致 alpha 过小
                PhaseFixOptions pfx;
                pfx.eig_floor = 1e-3; 
                pfx.inv_tol   = 1e-12;
                
                fixPhaseHessian(sys.phases[j], pfx);

                Ms[j] = sys.phases[j].M;
                mus[j] = sys.phases[j].mu;
                nPhases[j] = sys.phases[j].state.moleNumbers;
            }

            assembleGlobalSystem(
                sys.temperature,
                Ms, mus, nPhases,
                sys.elementMatrix,
                sys.feedMoles,
                Acoef, rhs
            );

            double lin_resid = 0.0;
            sol = solveGlobalLinearSystem(Acoef, rhs, &lin_resid);
            
            for(size_t i = E; i < E + F ; i++){
                std::cout<<"∆β["<<i-E<<"]="<<sol[i]<<"  ";
            }
            std::cout<<std::endl;

            // 抑制虚假的高残差警告（相对误差可能很小）
            // 只有当残差真的非常大时才打印
            // (feed - current) 产生的 RHS 项可能很大，导致 lin_resid 绝对值大，这是正常的
            if (lin_resid > 1e-4 * (1.0 + std::abs(sol[0]))) { // 简单的相对检查
                 // std::cerr << "[RandFlash] Linear resid: " << lin_resid << "\n";
            }

            if (sol.size() != E + F) throw std::runtime_error("Linear solve size mismatch");

            Lambda.assign(sol.begin(), sol.begin() + E);
            for(size_t j=0; j<F; ++j) deltaBeta[j] = sol[E + j];

            backSubstituteDeltas(
                sys.temperature,
                sys.elementMatrix,
                Ms, mus, nPhases,
                Lambda, deltaBeta,
                dnPhases
            );

            std::vector<std::vector<double>> currentNs(F);
            for(size_t j=0; j<F; ++j) currentNs[j] = sys.phases[j].state.moleNumbers;

            applyUpdate(sys.temperature, mus, dnPhases, currentNs);

            for(size_t j=0; j<F; ++j) sys.phases[j].state.moleNumbers = currentNs[j];

            auto conv = checkConvergence(
                iter + 1, tol, sys.temperature,
                mus, sys.elementMatrix, currentNs, sys.feedMoles
            );

            if (conv.converged) {
                result.success = true;
                result.iterations = iter + 1;
                result.mu_infinity_norm = conv.max_mu_diff;
                result.elem_residual_inf = conv.elem_error;
                
                result.n_phase = currentNs;
                result.beta.resize(F);
                for(size_t j=0; j<F; ++j) 
                    result.beta[j] = std::accumulate(currentNs[j].begin(), currentNs[j].end(), 0.0);
                
                return result;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[RandFlash] Exception: " << e.what() << std::endl;
        result.success = false;
        return result;
    }

    result.success = false;
    return result;
}