// rand_init.cpp
// 负责：两相/三相通用初始化（Least-Squares + Strict Element Scaling），生成可行初值（n,x,beta）。

#include "rand_flash.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace randflash {

// ---- moved here: two-phase initialization ----
InitResult RandFlash::initializeTwoPhase(
  const SystemContext& sys,
  const std::vector<double>& vaporGuess,
  const std::vector<double>& liquidGuess,
  double margin) const
{
  const size_t C = sys.feedMoles.size();
  auto sum = [](const std::vector<double>& v){ return std::accumulate(v.begin(), v.end(), 0.0); };
  const double Ftot = sum(sys.feedMoles);
  if (Ftot <= 0.0) throw std::invalid_argument("Feed total must be positive");

  // 本地辅助：归一化
  auto normalized = [](std::vector<double> v){
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    if (!(s > 0.0)) throw std::invalid_argument("Composition sum must be positive");
    for (double &xi : v) if (xi < 0.0) xi = 0.0;
    double s2 = std::accumulate(v.begin(), v.end(), 0.0);
    for (double &x : v) x /= (s2 > 0.0 ? s2 : 1.0);
    return v;
  };

  // 进料分率 z
  std::vector<double> z(C);
  for (size_t i=0;i<C;++i) z[i] = sys.feedMoles[i] / Ftot;

  // 准备通用输出结构：固定为2相
  InitResult out;
  out.n_phases.resize(2);       // index 0: Vapor, 1: Liquid
  out.beta.resize(2);
  out.compositions.resize(2);

  // 使用本地变量暂存计算结果，最后填入 out
  std::vector<double> y(C, 0.0), x(C, 0.0); // y=vapor, x=liquid
  double betaV = 0.0;

  const double eps = 1e-14;
  const double Fmin = eps, Fmax = Ftot * (1.0 - margin);

  // CASE A：只给了气相组成 y
  if (!vaporGuess.empty() && liquidGuess.empty()) {
    y = normalized(vaporGuess);

    double ub = Fmax;
    for (size_t i=0;i<C;++i) {
      if (y[i] > eps) ub = std::min(ub, sys.feedMoles[i] / y[i]);
    }
    ub = std::max(ub, Fmin);
    betaV = std::min(0.5*Ftot, ub*(1.0 - margin));

    const double Ltot = Ftot - betaV;
    if (Ltot <= Ftot*margin) betaV = Ftot*(1.0 - margin);
    
    for (size_t i=0;i<C;++i) {
      x[i] = (sys.feedMoles[i] - betaV*y[i]) / (Ftot - betaV);
      if (x[i] < 0.0) x[i] = 0.0;
    }
    x = normalized(x);
  }
  // CASE B：只给了液相组成 x
  else if (vaporGuess.empty() && !liquidGuess.empty()) {
    x = normalized(liquidGuess);

    double lb = Fmin;
    for (size_t i=0;i<C;++i) {
      if (x[i] > eps) lb = std::max(lb, Ftot - sys.feedMoles[i]/x[i]);
    }
    lb = std::max(lb, Fmin);
    betaV = std::max(0.5*Ftot, lb*(1.0 + margin));
    if (betaV >= Fmax) betaV = 0.5*(lb + Fmax);

    for (size_t i=0;i<C;++i) {
      y[i] = (sys.feedMoles[i] - (Ftot - betaV)*x[i]) / betaV;
      if (y[i] < 0.0) y[i] = 0.0;
    }
    y = normalized(y);
  }
  // CASE C：同时给了 y 与 x
  else if (!vaporGuess.empty() && !liquidGuess.empty()) {
    y = normalized(vaporGuess);
    x = normalized(liquidGuess);

    double num = 0.0, den = 0.0;
    for (size_t i=0;i<C;++i) {
      const double d = y[i] - x[i];
      num += d * (z[i] - x[i]);
      den += d * d;
    }
    betaV = (den > 0.0) ? Ftot * (num / den) : 0.5*Ftot;

    double ub = Fmax, lb = Fmin;
    for (size_t i=0;i<C;++i) {
      if (y[i] > eps) ub = std::min(ub, sys.feedMoles[i] / y[i]);
      if (x[i] > eps) lb = std::max(lb, Ftot - sys.feedMoles[i]/x[i]);
    }
    betaV = std::min(std::max(betaV, lb*(1.0 + margin)), ub*(1.0 - margin));
  }
  // CASE D：都没给 -> Wilson K
  else {
    std::vector<double> K(C, 0.0);
    // 这里要注意：如果 initializeTwoPhase 被用于非 0/1 相的情况，逻辑需调整
    // 但目前它是专门给 solveTwoPhase 用的
    thermo_.wilsonK(sys.temperature, sys.pressure, K);

    auto rr = [&](double b) {
      double s = 0.0;
      for (size_t i = 0; i < C; ++i) {
        double denom = 1.0 + b * (K[i] - 1.0);
        if (denom < 1e-12) denom = 1e-12;
        s += z[i] * (K[i] - 1.0) / denom;
      }
      return s;
    };

    // 简单的 Rachford-Rice 找 beta (0~1)
    double b_lo = margin, b_hi = 1.0 - margin;
    double f_lo = rr(b_lo), f_hi = rr(b_hi);
    double b_frac = 0.5;

    if (f_lo * f_hi < 0.0) {
      for (int it = 0; it < 40; ++it) {
        b_frac = 0.5 * (b_lo + b_hi);
        double f = rr(b_frac);
        if (std::fabs(f) < 1e-12) break;
        if (f * f_lo > 0.0) { b_lo = b_frac; f_lo = f; } 
        else { b_hi = b_frac; f_hi = f; }
      }
    }

    betaV = b_frac * Ftot;

    for (size_t i = 0; i < C; ++i) {
      double denom = 1.0 + b_frac * (K[i] - 1.0);
      if (denom < 1e-12) denom = 1e-12;
      x[i] = z[i] / denom;
      y[i] = K[i] * x[i];
    }
    x = normalized(x);
    y = normalized(y);
  }

  // 填回通用结构
  double betaL = Ftot - betaV;
  
  out.beta[0] = betaV;
  out.beta[1] = betaL;
  
  out.compositions[0] = y;
  out.compositions[1] = x;

  out.n_phases[0].resize(C);
  out.n_phases[1].resize(C);
  for (size_t i=0; i<C; ++i) {
    out.n_phases[0][i] = betaV * y[i];
    out.n_phases[1][i] = betaL * x[i];
  }

  return out;
}


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

// ========== Reactive System Initialization ==========

InitResult RandFlash::initializeReactiveSinglePhase(
    const SystemContext& sys,
    const std::vector<double>& elementMoles,
    int phaseFlag) const
{
    const size_t C = sys.feedMoles.size();
    const size_t E = sys.elementMatrix.size();

    if (elementMoles.size() != E) {
        throw std::invalid_argument(
            "Element moles size (" + std::to_string(elementMoles.size()) +
            ") does not match element matrix rows (" + std::to_string(E) + ")");
    }

    // Solve least-squares problem: min ||n||^2 s.t. A·n = b, n ≥ 0
    // Using simple approach: start with uniform distribution, then project to satisfy A·n = b

    std::vector<double> n(C);

    if (E == C) {
        // Square system: try direct solve A·n = b
        // For simplicity, use pseudo-inverse approach via least squares
        // This is a simplified implementation - could use Eigen for better numerics

        // Start with uniform guess
        double total = std::accumulate(elementMoles.begin(), elementMoles.end(), 0.0);
        for (size_t i = 0; i < C; ++i) {
            n[i] = total / static_cast<double>(C);
        }
    } else {
        // Underdetermined system (E < C): use minimum norm solution
        // Start with a feasible guess based on element ratios

        double total = std::accumulate(elementMoles.begin(), elementMoles.end(), 0.0);
        for (size_t i = 0; i < C; ++i) {
            n[i] = total / static_cast<double>(C);
        }
    }

    // Iterative projection to satisfy element conservation
    // This ensures A·n = b while keeping n ≥ 0
    const int max_iter = 100;
    const double tol = 1e-12;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Compute current element moles: b_computed = A·n
        std::vector<double> b_computed(E, 0.0);
        for (size_t e = 0; e < E; ++e) {
            for (size_t i = 0; i < C; ++i) {
                b_computed[e] += sys.elementMatrix[e][i] * n[i];
            }
        }

        // Check convergence
        double max_error = 0.0;
        for (size_t e = 0; e < E; ++e) {
            double error = std::abs(b_computed[e] - elementMoles[e]);
            max_error = std::max(max_error, error);
        }

        if (max_error < tol) break;

        // Correction step: adjust n to reduce element error
        // For each element constraint, distribute the error proportionally
        for (size_t e = 0; e < E; ++e) {
            double error = elementMoles[e] - b_computed[e];

            // Find species containing this element
            double sum_A = 0.0;
            for (size_t i = 0; i < C; ++i) {
                if (sys.elementMatrix[e][i] > 1e-10) {
                    sum_A += sys.elementMatrix[e][i];
                }
            }

            if (sum_A > 1e-10) {
                // Distribute error proportionally to element coefficients
                for (size_t i = 0; i < C; ++i) {
                    if (sys.elementMatrix[e][i] > 1e-10) {
                        n[i] += error * sys.elementMatrix[e][i] / sum_A;
                        n[i] = std::max(n[i], 1e-20);  // Keep positive
                    }
                }
            }
        }
    }

    // Ensure all moles are positive
    for (double& ni : n) {
        ni = std::max(ni, 1e-20);
    }

    // Compute total moles and composition
    double n_total = std::accumulate(n.begin(), n.end(), 0.0);
    std::vector<double> x(C);
    for (size_t i = 0; i < C; ++i) {
        x[i] = n[i] / n_total;
    }

    // Build result
    InitResult result;
    result.n_phases.resize(1);
    result.beta.resize(1);
    result.compositions.resize(1);

    result.n_phases[0] = n;
    result.beta[0] = n_total;
    result.compositions[0] = x;

    return result;
}

InitResult RandFlash::initializeReactiveMultiPhase(
    const SystemContext& sys,
    int numPhases,
    const std::vector<int>& phaseFlags,
    const std::vector<double>& elementMoles) const
{
    const size_t C = sys.feedMoles.size();
    const size_t E = sys.elementMatrix.size();

    if (numPhases < 1) {
        throw std::invalid_argument("Number of phases must be at least 1");
    }

    if (phaseFlags.size() != static_cast<size_t>(numPhases)) {
        throw std::invalid_argument(
            "Phase flags size (" + std::to_string(phaseFlags.size()) +
            ") does not match number of phases (" + std::to_string(numPhases) + ")");
    }

    if (elementMoles.size() != E) {
        throw std::invalid_argument(
            "Element moles size does not match element matrix rows");
    }

    // Strategy: distribute elements among phases based on phase type
    // For simplicity, start with equal distribution, then use Wilson K for vapor/liquid split

    InitResult result;
    result.n_phases.resize(numPhases);
    result.beta.resize(numPhases);
    result.compositions.resize(numPhases);

    // Total element moles
    double total_elements = std::accumulate(elementMoles.begin(), elementMoles.end(), 0.0);

    // Initial guess: distribute elements equally among phases
    std::vector<double> phase_element_fraction(numPhases, 1.0 / numPhases);

    // For each phase, compute initial composition
    for (int j = 0; j < numPhases; ++j) {
        // Allocate fraction of elements to this phase
        std::vector<double> phase_element_moles(E);
        for (size_t e = 0; e < E; ++e) {
            phase_element_moles[e] = elementMoles[e] * phase_element_fraction[j];
        }

        // Create temporary system context for single-phase initialization
        SystemContext phase_sys = sys;

        // Initialize this phase composition from its element allocation
        auto phase_init = initializeReactiveSinglePhase(phase_sys, phase_element_moles, phaseFlags[j]);

        result.n_phases[j] = phase_init.n_phases[0];
        result.beta[j] = phase_init.beta[0];
        result.compositions[j] = phase_init.compositions[0];
    }

    // Verify total element conservation
    std::vector<double> total_computed(E, 0.0);
    for (int j = 0; j < numPhases; ++j) {
        for (size_t e = 0; e < E; ++e) {
            for (size_t i = 0; i < C; ++i) {
                total_computed[e] += sys.elementMatrix[e][i] * result.n_phases[j][i];
            }
        }
    }

    // Scale to match exact element conservation
    for (size_t e = 0; e < E; ++e) {
        if (total_computed[e] > 1e-20) {
            double scale = elementMoles[e] / total_computed[e];
            // Apply scaling to all phases proportionally
            for (int j = 0; j < numPhases; ++j) {
                for (size_t i = 0; i < C; ++i) {
                    if (sys.elementMatrix[e][i] > 1e-10) {
                        result.n_phases[j][i] *= scale;
                    }
                }
            }
        }
    }

    // Recompute beta and compositions after scaling
    for (int j = 0; j < numPhases; ++j) {
        result.beta[j] = std::accumulate(result.n_phases[j].begin(), result.n_phases[j].end(), 0.0);
        for (size_t i = 0; i < C; ++i) {
            result.compositions[j][i] = result.n_phases[j][i] / result.beta[j];
        }
    }

    return result;
}


} // namespace randflash
