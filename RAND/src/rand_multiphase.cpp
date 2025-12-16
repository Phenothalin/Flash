#include "rand_flash.hpp"
#include "phase_stability.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
using namespace randflash;
using namespace ls;

namespace {

static int find_water_index(const thermo::IThermoBackend& thermo, size_t C) {
    const auto names = thermo.getComponentNames();
    for (size_t i = 0; i < std::min(names.size(), C); ++i) {
        std::string n = names[i];
        std::transform(n.begin(), n.end(), n.begin(), ::toupper);
        if (n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

static void reorder_flash_result(
    MultiFlashResult& res,
    const std::vector<int>& phaseFlags,
    int vaporFlag,
    int liquidFlag,
    int water_idx)
{
    const size_t F = res.n_phase.size();
    if (F <= 1 || phaseFlags.size() != F) return;

    // Build (idx, isVap, beta, x_water)
    struct Key {
        size_t idx;
        bool isVap;
        double beta;
        double xw;
    };
    std::vector<Key> keys;
    keys.reserve(F);

    for (size_t j = 0; j < F; ++j) {
        double beta = 0.0;
        for (double v : res.n_phase[j]) beta += v;
        double xw = 0.0;
        if (water_idx >= 0 && beta > 0.0 && static_cast<size_t>(water_idx) < res.n_phase[j].size()) {
            xw = res.n_phase[j][static_cast<size_t>(water_idx)] / beta;
        }
        keys.push_back(Key{j, phaseFlags[j] == vaporFlag, beta, xw});
    }

    // Primary sort: vapor first
    // Secondary: for liquids, if water exists and there are >=2 liquids, order by water fraction ascending
    // Otherwise: by beta descending
    const bool hasWater = (water_idx >= 0);
    const int nLiq = static_cast<int>(std::count_if(keys.begin(), keys.end(), [](const Key& k){ return !k.isVap; }));

    std::sort(keys.begin(), keys.end(), [&](const Key& a, const Key& b){
        if (a.isVap != b.isVap) return a.isVap > b.isVap;
        if (!a.isVap && !b.isVap && hasWater && nLiq >= 2) {
            // oil first (low xw), water-rich later
            if (std::abs(a.xw - b.xw) > 1e-8) return a.xw < b.xw;
        }
        return a.beta > b.beta;
    });

    // Apply permutation
    std::vector<std::vector<double>> n_new;
    std::vector<double> beta_new;
    n_new.reserve(F);
    beta_new.reserve(F);

    for (const auto& k : keys) {
        n_new.push_back(res.n_phase[k.idx]);
        beta_new.push_back(res.beta.empty() ? k.beta : res.beta[k.idx]);
    }
    res.n_phase.swap(n_new);
    res.beta.swap(beta_new);
}

} // anonymous namespace

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

// ----------------------------------------------------------------------------------
// 通用：由“相组成猜测”生成满足严格守恒的初始 n（最小二乘 beta + 行缩放）
// - phaseCompositions: F x C（可不严格归一，内部会归一）
// - phaseFlags: F (可为空；为空则默认 0 为气相，其余为液相)
// ----------------------------------------------------------------------------------
InitResult RandFlash::initializeFromCompositions(
    const SystemContext& sys,
    int nPhases,
    const std::vector<std::vector<double>>& phaseCompositions,
    const std::vector<int>& phaseFlags,
    double min_phase_moles_ratio) const
{
    const size_t C = sys.feedMoles.size();
    if (nPhases <= 0) throw std::invalid_argument("initializeFromCompositions: nPhases <= 0");
    if ((int)phaseCompositions.size() != nPhases) {
        throw std::invalid_argument("initializeFromCompositions: phaseCompositions.size() != nPhases");
    }

    const double total_feed = std::accumulate(sys.feedMoles.begin(), sys.feedMoles.end(), 0.0);

    // 1) 归一化组成
    std::vector<std::vector<double>> X_norm(nPhases, std::vector<double>(C, 0.0));
    for (int j = 0; j < nPhases; ++j) {
        double s = (j < (int)phaseCompositions.size())
                 ? std::accumulate(phaseCompositions[j].begin(), phaseCompositions[j].end(), 0.0)
                 : 0.0;
        if (s < 1e-30) s = 1.0;
        for (size_t i = 0; i < C; ++i) {
            double v = (i < phaseCompositions[j].size()) ? phaseCompositions[j][i] / s : 0.0;
            if (!std::isfinite(v) || v < 0.0) v = 0.0;
            X_norm[j][i] = v;
        }
        // 若全 0，退化为均匀
        double s2 = std::accumulate(X_norm[j].begin(), X_norm[j].end(), 0.0);
        if (s2 <= 0.0) {
            const double uni = 1.0 / static_cast<double>(C);
            std::fill(X_norm[j].begin(), X_norm[j].end(), uni);
        } else {
            for (double& v : X_norm[j]) v /= s2;
        }
    }

    // 2) 正规方程 A * beta = b
    std::vector<double> A_ls(nPhases * nPhases, 0.0);
    std::vector<double> b_ls(nPhases, 0.0);
    for (int j = 0; j < nPhases; ++j) {
        for (size_t i = 0; i < C; ++i) {
            b_ls[j] += X_norm[j][i] * sys.feedMoles[i];
        }
        for (int k = 0; k < nPhases; ++k) {
            double dot = 0.0;
            for (size_t i = 0; i < C; ++i) dot += X_norm[j][i] * X_norm[k][i];
            A_ls[j * nPhases + k] = dot;
        }
    }

    double dummy_res = 0.0;
    std::vector<double> betas = linearSolver_.solveDense(nPhases, A_ls, b_ls, &dummy_res);

    // 3) beta 下限保护
    const double min_beta = std::max(1e-12, min_phase_moles_ratio * total_feed);
    for (double& b : betas) {
        if (!std::isfinite(b) || b < min_beta) b = min_beta;
    }

    // 4) 初步 n = beta * x，然后按组分行缩放保证严格守恒
    std::vector<std::vector<double>> n_temp(nPhases, std::vector<double>(C, 0.0));
    std::vector<double> n_sum(C, 0.0);
    for (int j = 0; j < nPhases; ++j) {
        for (size_t i = 0; i < C; ++i) {
            n_temp[j][i] = betas[j] * X_norm[j][i];
            n_sum[i] += n_temp[j][i];
        }
    }
    for (size_t i = 0; i < C; ++i) {
        const double target = sys.feedMoles[i];
        const double current = n_sum[i];
        if (current > 1e-30) {
            const double scale = target / current;
            for (int j = 0; j < nPhases; ++j) n_temp[j][i] *= scale;
        } else {
            // 所有相对该组分都为 0，但进料不为 0：平均分配
            for (int j = 0; j < nPhases; ++j) n_temp[j][i] = target / static_cast<double>(nPhases);
        }
    }

    // 5) 输出
    InitResult out;
    out.n_phases = n_temp;
    out.beta.resize(nPhases, 0.0);
    out.compositions.resize(nPhases, std::vector<double>(C, 0.0));
    for (int j = 0; j < nPhases; ++j) {
        out.beta[j] = std::accumulate(out.n_phases[j].begin(), out.n_phases[j].end(), 0.0);
        out.compositions[j] = out.n_phases[j];
        if (out.beta[j] > 1e-30) {
            const double invB = 1.0 / out.beta[j];
            for (double& v : out.compositions[j]) v *= invB;
        } else {
            const double uni = 1.0 / static_cast<double>(C);
            std::fill(out.compositions[j].begin(), out.compositions[j].end(), uni);
        }
    }

    (void)phaseFlags; // phaseFlags 由调用处决定如何赋给 PhaseState
    return out;
}

// ----------------------------------------------------------------------------------
// 通用三相初始化：
// 1) 若检测到 H2O 且含量显著，优先使用顺序拆分（水相/油相/气相）
// 2) 否则：基于 Wilson-K 构造 Vapor / Heavy-Liquid / Intermediate-Liquid 三个组成猜测
// 再用 initializeFromCompositions 生成严格守恒的 n。
// ----------------------------------------------------------------------------------
InitResult RandFlash::initializeThreePhaseGeneric(const SystemContext& sys, double margin) const
{
    const size_t C = sys.feedMoles.size();
    if (C == 0) return {};

    // 归一化 z
    const double Ftot = std::accumulate(sys.feedMoles.begin(), sys.feedMoles.end(), 0.0);
    std::vector<double> z = sys.feedMoles;
    if (Ftot > 1e-30) {
        for (double& v : z) v /= Ftot;
    } else {
        const double uni = 1.0 / static_cast<double>(C);
        std::fill(z.begin(), z.end(), uni);
    }

    // 识别水
    auto names = thermo_.getComponentNames();
    int water_idx = -1;
    for (size_t i = 0; i < std::min(names.size(), C); ++i) {
        std::string n = names[i];
        std::transform(n.begin(), n.end(), n.begin(), ::toupper);
        if (n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
            water_idx = static_cast<int>(i);
            break;
        }
    }
    const double z_water = (water_idx >= 0) ? z[water_idx] : 0.0;
    if (water_idx >= 0 && z_water > 1e-4) {
        // 直接复用你已经验证过的含水初始化
        return initializeThreePhaseWater(sys, margin);
    }

    // Wilson K
    std::vector<double> K(C, 1.0);
    thermo_.wilsonK(sys.temperature, sys.pressure, K);
    for (double& Ki : K) Ki = std::clamp(Ki, 1e-8, 1e8);

    auto normalize = [&](std::vector<double> x) {
        for (double& v : x) v = std::max(v, margin);
        double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s <= 0.0) {
            const double uni = 1.0 / static_cast<double>(C);
            std::fill(x.begin(), x.end(), uni);
            return x;
        }
        for (double& v : x) v /= s;
        return x;
    };

    // 三个组成猜测
    std::vector<double> xV(C), xLh(C), xLi(C);
    for (size_t i = 0; i < C; ++i) {
        xV[i]  = z[i] * K[i];                 // vapor-like
        xLh[i] = z[i] / K[i];                 // heavy liquid-like
        xLi[i] = z[i] / std::sqrt(K[i]);      // intermediate liquid-like
    }
    xV  = normalize(xV);
    xLh = normalize(xLh);
    xLi = normalize(xLi);

    std::vector<std::vector<double>> comps = {xV, xLh, xLi};
    std::vector<int> flags = {thermo_.vaporPhaseFlag(), thermo_.liquidPhaseFlag(), thermo_.liquidPhaseFlag()};
    return initializeFromCompositions(sys, 3, comps, flags);
}

  MultiFlashResult RandFlash::solveMultiPhaseCore(
  const FlashInput& input,
  int nPhases,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<std::vector<double>>& initialPhaseCompositions,
  const std::vector<int>& initialPhaseFlags,
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
        std::vector<int> flags = initialPhaseFlags;
        if ((int)flags.size() != 2) {
            flags = {thermo_.vaporPhaseFlag(), thermo_.liquidPhaseFlag()};
        }
        sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], flags[0]};
        sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], flags[1]};
        initialized = true;
    } 
    // 情况 B: 三相且无初值 -> 调用通用三相初始化
    else if (nPhases == 3 && initialPhaseCompositions.empty()) {
        auto init = initializeThreePhaseGeneric(sys);
        std::vector<int> flags = initialPhaseFlags;
        if ((int)flags.size() != 3) {
            flags = {thermo_.vaporPhaseFlag(), thermo_.liquidPhaseFlag(), thermo_.liquidPhaseFlag()};
        }
        sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], flags[0]};
        sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], flags[1]};
        sys.phases[2].state = {input.temperature, input.pressure, init.n_phases[2], flags[2]};
        initialized = true;
    }
    // 情况 C: 有用户提供的初值 -> 最小二乘投影 + 强制守恒缩放
    else if (!initialPhaseCompositions.empty()) {
        if ((int)initialPhaseCompositions.size() != nPhases) {
            std::cerr << "[Warning] Initial compositions count != nPhases. Using default logic.\n";
        } else {
            std::vector<int> flags = initialPhaseFlags;
            if ((int)flags.size() != nPhases) {
                flags.assign(nPhases, thermo_.liquidPhaseFlag());
                flags[0] = thermo_.vaporPhaseFlag();
            }

            auto init = initializeFromCompositions(sys, nPhases, initialPhaseCompositions, flags);
            for (int j = 0; j < nPhases; ++j) {
                sys.phases[j].state = {input.temperature, input.pressure, init.n_phases[j], flags[j]};
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
        std::vector<int> flags = initialPhaseFlags;
        if ((int)flags.size() != nPhases) {
            flags.assign(nPhases, thermo_.liquidPhaseFlag());
            flags[0] = thermo_.vaporPhaseFlag();
        }
        for(int j=0; j<nPhases; ++j) {
            int flag = flags[j];
            sys.phases[j].state = {input.temperature, input.pressure, n, flag};
        }
    }

    // 3. 核心求解
    auto res = solveGeneral(sys, maxIter, tol);

    // 4. 结果相序统一：Vapor -> (Oil-like liquid) -> (Water-like liquid)
    // 说明：算法内部相序并不要求固定，但测试与输出更希望保持稳定顺序，
    // 否则会造成“油/水相对调”的观感，并可能影响外部调用。
    if (res.success && res.n_phase.size() == sys.phases.size()) {
        std::vector<int> flags;
        flags.reserve(sys.phases.size());
        for (const auto& ph : sys.phases) flags.push_back(ph.state.phaseFlag);
        const int vap = thermo_.vaporPhaseFlag();
        const int liq = thermo_.liquidPhaseFlag();
        const int water_idx = find_water_index(thermo_, input.feedMoles.size());
        reorder_flash_result(res, flags, vap, liq, water_idx);
    }

    return res;
}

// ----------------------------------------------------------------------------------
// 通用接口：自动相稳定性分析 -> 自动相数选择 -> 自动初始化 -> 相分裂
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solve(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions,
    int maxIter,
    double tol)
{
    // 若用户给了初值，则按给定相数直接求解（相标志保持默认：0 为气相，其余为液相）
    if (!initialPhaseCompositions.empty()) {
        const int F = static_cast<int>(initialPhaseCompositions.size());
        return solveMultiPhaseCore(input, F, elementMatrix, initialPhaseCompositions, {}, maxIter, tol);
    }

    // 1) 相稳定性分析
    phase_stability::PhaseStabilityAnalyzer analyzer(thermo_);
    phase_stability::StabilityOptions opt;
    opt.verbose = false;

    auto stab = analyzer.analyze(input.temperature, input.pressure, input.feedMoles, opt);

    // 2) 单相稳定：直接返回（不进入 Rand 迭代）
    if (stab.stable) {
        MultiFlashResult result;
        result.success = true;
        result.iterations = 0;
        result.pressure = input.pressure;
        result.temperature = input.temperature;
        result.feedComposition = input.feedMoles;

        const double tot = std::accumulate(input.feedMoles.begin(), input.feedMoles.end(), 0.0);
        result.n_phase = {input.feedMoles};
        result.beta = {tot};
        result.mu_infinity_norm = 0.0;
        result.elem_residual_inf = 0.0;
        return result;
    }

    // 3) 根据稳定性结果选择相数与相类型。
    // 关键约束：对“典型 VLE”体系（一个气相+一个液相）
    // 稳定性分析可能同时给出 vap 和 liq 的 incipient（数值噪声 / seed 多样性导致），
    // 但这不应被直接解释为 3 相。3 相仅在“需要两套液相”（LLE/LLV）时启用。

    auto normalize = [&](const std::vector<double>& v) {
        std::vector<double> x = v;
        for (double& xi : x) xi = std::max(xi, 1e-14);
        double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s <= 0.0) {
            const double uni = 1.0 / static_cast<double>(x.size());
            std::fill(x.begin(), x.end(), uni);
        } else {
            for (double& xi : x) xi /= s;
        }
        return x;
    };

    const int vap = thermo_.vaporPhaseFlag();
    const int liq = thermo_.liquidPhaseFlag();

    // Wilson K (用于在 stability 未给出某类 incipient 时构造更“有偏”的相组成初猜)
    std::vector<double> K(input.feedMoles.size(), 1.0);
    try {
        thermo_.wilsonK(input.temperature, input.pressure, K);
        for (double& Ki : K) Ki = std::clamp(Ki, 1e-12, 1e12);
    } catch (...) {
        // fallback: keep K=1
    }

    auto build_wilson_guess = [&](bool vapor_like) {
        std::vector<double> x = input.feedMoles;
        const double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s > 0.0) {
            for (double& xi : x) xi /= s;
        }
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = vapor_like ? (x[i] * K[i]) : (x[i] / K[i]);
        }
        return normalize(x);
    };

    // 探测是否含水（用于 LLV：水/油两液相）
    int water_idx = -1;
    {
        const auto names = thermo_.getComponentNames();
        for (size_t i = 0; i < std::min(names.size(), input.feedMoles.size()); ++i) {
            std::string n = names[i];
            std::transform(n.begin(), n.end(), n.begin(), ::toupper);
            if (n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
                water_idx = static_cast<int>(i);
                break;
            }
        }
    }
    const double ztot = std::accumulate(input.feedMoles.begin(), input.feedMoles.end(), 0.0);
    const double z_water = (water_idx >= 0 && ztot > 0.0) ? (input.feedMoles[water_idx] / ztot) : 0.0;
    const bool has_water = (water_idx >= 0 && z_water > 1e-6);

    // 收集 incipient：最多 2 个
    std::vector<std::vector<double>> liquid_cands;
    std::vector<std::vector<double>> vapor_cands;
    for (const auto& inc : stab.incipient) {
        if (inc.phase_flag == liq) liquid_cands.push_back(inc.x);
        if (inc.phase_flag == vap) vapor_cands.push_back(inc.x);
    }

    // 默认：先做 2 相（最不容易引入回归）
    // VLE: {VAP, LIQ}
    // LLE: {LIQ, LIQ}
    // LLV: {VAP, LIQ, LIQ}

    // --- 2 相初猜 ---
    std::vector<std::vector<double>> comps2;
    std::vector<int> flags2;
    comps2.reserve(2);
    flags2.reserve(2);

    // vapor guess
    std::vector<double> xV = (!vapor_cands.empty()) ? vapor_cands.front() : build_wilson_guess(true);
    // liquid guess
    std::vector<double> xL = (!liquid_cands.empty()) ? liquid_cands.front() : build_wilson_guess(false);

    // 若 reference 是 vapor 或 liquid，则更倾向于把 z 放在对应相上（更稳健）
    if (stab.reference_phase_flag == vap) xV = normalize(input.feedMoles);
    if (stab.reference_phase_flag == liq) xL = normalize(input.feedMoles);

    // 判断是否更像 LLE（两套 liquid），仅当 stability 明确给出 >=2 个 liquid candidate
    const bool want_LLE = (liquid_cands.size() >= 2) && (vapor_cands.empty()) && !has_water;

    if (want_LLE) {
        comps2 = { liquid_cands[0], liquid_cands[1] };
        flags2 = { liq, liq };
    } else {
        comps2 = { xV, xL };
        flags2 = { vap, liq };
    }

    // --- 3 相判别（仅当确实需要两套 liquid 时才启用）---
    // 1) 水体系：优先尝试 LLV
    // 2) 非水体系：只有在 liquid candidate >= 2 时才尝试 LLV
    const bool want_LLV = has_water || (liquid_cands.size() >= 2);

    if (want_LLV && !want_LLE) {
        // 3 相初值：优先来自稳定性结果；不够则用通用 3 相初始化补齐
        SystemContext sys0;
        sys0.temperature = input.temperature;
        sys0.pressure = input.pressure;
        sys0.feedMoles = input.feedMoles;
        sys0.elementMatrix = elementMatrix;

        std::vector<std::vector<double>> comps3;
        std::vector<int> flags3 = {vap, liq, liq};

        // vapor
        if (!vapor_cands.empty()) comps3.push_back(vapor_cands.front());
        else comps3.push_back(build_wilson_guess(true));

        // two liquids
        if (liquid_cands.size() >= 2) {
            comps3.push_back(liquid_cands[0]);
            comps3.push_back(liquid_cands[1]);
        } else {
            // 用 Generic 初始化构造两套液相种子（尤其对含水体系更稳健）
            auto init3 = initializeThreePhaseGeneric(sys0);
            // init3.compositions 是 {V, L1, L2}
            comps3 = init3.compositions;
        }

        auto res3 = solveMultiPhaseCore(input, 3, elementMatrix, comps3, flags3, maxIter, tol);
        if (res3.success) return res3;
        // 3 相失败时回退 2 相
    }

    return solveMultiPhaseCore(input, 2, elementMatrix, comps2, flags2, maxIter, tol);
}

MultiFlashResult RandFlash::solveMultiPhase(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions,
    const std::vector<int>& initialPhaseFlags,
    int maxIter,
    double tol)
{
  // 兼容旧接口：若用户显式给了相组成，则按给定相数求解；否则走“自动相稳定性分析”。
  if (!initialPhaseCompositions.empty()) {
    const int F = static_cast<int>(initialPhaseCompositions.size());
    return solveMultiPhaseCore(
        input,
        F,
        elementMatrix,
        initialPhaseCompositions,
        initialPhaseFlags,
        maxIter,
        tol);
  }

  return solve(input, elementMatrix, {}, maxIter, tol);
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
                
                PhaseFixOptions pfx;
                // 说明：这里不应把 eig_floor 设得过大（例如 1e-3），否则会过度“钝化”真实的曲率信息，
                // 导致步长被压得很小、迭代次数增加，甚至破坏原本可收敛的两相工况。
                // 维持与两相版本一致的默认量级，并让 line-search 决定实际步长。
                pfx.eig_floor = 1e-10;
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