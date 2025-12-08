// RAND/src/rand_flash.cpp
#include "rand_flash.hpp"
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace randflash;
using namespace ls;


// === 在切空间做 SPD 修正，并保持 m x = 1 (使用 Woodbury 公式避免病态求逆) ===
// 输入: solver, x, m_in (CxC)
// 输出: PhaseFixResult (包含修正后的 M 和 m)
static PhaseFixResult fix_phase_hessian_one_phase(
  ls::LinearSolverInterface& solver,
  const std::vector<double>& x,                 
  const std::vector<std::vector<double>>& m_in, 
  const PhaseFixOptions& opt = {})
{
  const int C = static_cast<int>(x.size());
  PhaseFixResult out;
  
  // 基础校验
  if (C == 0) {
    out.applied = false; return out;
  }

  // 1. 构建 Scaled Hessian Q 以及 Scaling 因子 D^{1/2}
  //    m = D^{-1/2} Q D^{-1/2}, 其中 D = diag(x)
  //    所以 Q_{ij} = sqrt(x_i) * m_{ij} * sqrt(x_j)
  //    Q 是良态的 (特征值通常在 1 附近)
  
  std::vector<double> sqrt_x(C);
  for(int i=0; i<C; ++i) sqrt_x[i] = sqrt(x[i]);

  std::vector<double> Q_flat(C * C);
  for(int i=0; i<C; ++i) {
      for(int j=0; j<C; ++j) {
          Q_flat[i*C + j] = sqrt_x[i] * m_in[i][j] * sqrt_x[j];
      }
  }

  // 2. 切空间投影与特征值分析 (基于 Q)
  //    我们需要检查 Q 在切空间 (constrain: sum n_i = const) 上的正定性
  //    切空间基 B (C x C-1)
  //    实际上 Q 的特征向量 u 对应 m 的特征向量 v = D^{-1/2} u
  
  // 这里为了简单和保持逻辑一致，我们沿用论文的方法：
  // 直接对 Q 矩阵进行特征值分解。因为 Q 是对称矩阵。
  // 注意：m 矩阵必有一个特征值为 1 对应特征向量 x (因为 m*x = 1)
  // 对应地，Q * sqrt(x) = sqrt(x)，因为:
  // (sqrt(x) m sqrt(x)) * sqrt(x) = sqrt(x) * (m x) = sqrt(x) * 1
  // 所以 Q 也有一个特征值为 1，特征向量为 sqrt(x)。这个方向是“非切空间”方向。
  // 我们只需要检查除这个方向以外的最小特征值。

  // 为了精确，还是使用切空间投影法
  auto B = ls::tangentBasis(C); // C x (C-1)
  const int T = C - 1;

  // Q_tan = B^T * Q * B
  std::vector<double> Q_tan_flat(T * T, 0.0);
  for(int p=0; p<T; ++p) {
      for(int q=0; q<T; ++q) {
          double val = 0.0;
          for(int i=0; i<C; ++i) {
              double temp = 0.0;
              for(int j=0; j<C; ++j) {
                  // Q_ij * B_jq
                  temp += Q_flat[i*C+j] * B[j][q];
              }
              val += B[i][p] * temp;
          }
          Q_tan_flat[p*T + q] = val;
      }
  }

  std::vector<double> evals;
  std::vector<double> evecs_sub; // 切空间特征向量
  solver.eigenDecomposeSymmetric(T, Q_tan_flat, evals, evecs_sub);

  if (evals.empty()) {
      // C=1 的情况，没有切空间，直接返回
      // M = x * x^T (单组分)
      out.applied = false;
      out.m_fixed = m_in;
      out.M_fixed.assign(C, std::vector<double>(C));
      out.M_fixed[0][0] = x[0]*x[0];
      return out;
  }

  double lam_min = evals[0];
  out.lam_min_before = lam_min;

  // 3. 判断是否需要修正
  if (lam_min >= opt.eig_floor) {
      // 不需要修正，直接计算 M = m^{-1}
      // 使用 Scaling 方法求逆以保证精度: M = D^{1/2} Q^{-1} D^{1/2}
      
      std::vector<double> Q_inv = solver.invertSPD(C, Q_flat);
      
      std::vector<std::vector<double>> M_calc(C, std::vector<double>(C));
      for(int i=0; i<C; ++i) {
          for(int j=0; j<C; ++j) {
              M_calc[i][j] = sqrt_x[i] * Q_inv[i*C+j] * sqrt_x[j];
          }
      }
      out.applied = false;
      out.lam_min_after = lam_min;
      out.m_fixed = m_in;
      out.M_fixed = M_calc;
      return out;
  }

  // 4. 需要修正：构造修正项
  //    目标: m_new = m + k * v * v^T
  //    其中 k = (target - lam_min), v 是 m 在切空间对应的特征向量
  //    注意特征向量转换: u_tan (T维) -> u_full (C维, Q的特征向量) -> v (C维, m的修正向量)
  
  double target = std::abs(lam_min) + opt.eig_floor + opt.eps_shift; // 修正量
  // 修正后的特征值设为略大于0的一个值，比如 0.001
  // 这里 k 是加在特征值上的量
  
  // 4.1 恢复 Q 的最小特征向量 u
  // u = B * evec_sub
  std::vector<double> u(C, 0.0);
  for(int i=0; i<C; ++i) {
      for(int k=0; k<T; ++k) {
          u[i] += B[i][k] * evecs_sub[0*T + k]; // evecs_sub 第一列是最小特征向量
      }
  }
  
  // 4.2 计算 m 的修正向量 v
  // 论文定义修正项为 k * u * u^T 加在 Q 上
  // 对应加在 m 上的项是 k * (D^{-1/2} u) * (D^{-1/2} u)^T
  // 令 v_i = u_i / sqrt(x_i)
  std::vector<double> v(C);
  for(int i=0; i<C; ++i) {
      v[i] = u[i] / sqrt_x[i];
  }

  // 5. 计算未修正的 M_raw (使用 Scaling 技巧)
  //    M_raw = D^{1/2} Q^{-1} D^{1/2}
  //    虽然 Q 有负特征值，但绝对值通常远离0，直接求逆通常是安全的（只要不是奇异）
  //    如果不放心，可以对 Q 先做修正再求逆，但使用 Woodbury 公式更优雅
  std::vector<double> Q_inv = solver.invertSPD(C, Q_flat);
  
  std::vector<std::vector<double>> M_raw(C, std::vector<double>(C));
  for(int i=0; i<C; ++i) {
      for(int j=0; j<C; ++j) {
          M_raw[i][j] = sqrt_x[i] * Q_inv[i*C+j] * sqrt_x[j];
      }
  }

  // 6. 使用 Woodbury 公式更新 M
  //    (m + k v v^T)^{-1} = M - (k (M v) (M v)^T) / (1 + k v^T M v)
  
  double k_val = target; 

  // 计算 temp = M_raw * v
  std::vector<double> Mv(C, 0.0);
  for(int i=0; i<C; ++i) {
      for(int j=0; j<C; ++j) {
          Mv[i] += M_raw[i][j] * v[j];
      }
  }

  // 计算分母 denom = 1 + k * v^T * M * v
  double vMv = 0.0;
  for(int i=0; i<C; ++i) vMv += v[i] * Mv[i];
  
  double denom = 1.0 + k_val * vMv;

  // 更新得到 M_fixed
  std::vector<std::vector<double>> M_fix(C, std::vector<double>(C));
  for(int i=0; i<C; ++i) {
      for(int j=0; j<C; ++j) {
          M_fix[i][j] = M_raw[i][j] - (k_val / denom) * Mv[i] * Mv[j];
      }
  }

  // 7. 更新 m_fixed (仅用于调试或完整性，主要用的是 M)
  //    m_fixed = m_in + k * v * v^T
  std::vector<std::vector<double>> m_fix = m_in;
  for(int i=0; i<C; ++i) {
      for(int j=0; j<C; ++j) {
          m_fix[i][j] += k_val * v[i] * v[j];
      }
  }

  out.applied = true;
  out.lam_min_after = lam_min + k_val; // 近似值
  out.m_fixed = m_fix;
  out.M_fixed = M_fix;

  return out;
}

RandFlash::RandFlash(thermo::IThermoBackend& thermo,
  ls::LinearSolverInterface& linearSolver)
  : thermo_(thermo),
  linearSolver_(linearSolver)
{
  // 如果有额外初始化,可以放在这里
}

// ---- in rand_flash.cpp ----
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

// 1) 局部 Jacobian 构造：对应论文式 (4.11)-(4.14)
void RandFlash::updatePhaseChemistry(PhaseContext& phaseCtx)
{
  // 1) 复制一份状态以便安全裁剪 n_i（避免零 / 负数）
  thermo::PhaseState st = phaseCtx.state;
  for (double& ni : st.moleNumbers) {
    if (ni <= 1e-10) ni = 1e-10;
  }

  const size_t C = st.moleNumbers.size();
  const double RT = R_CONST * st.Temperature;

  // 2) 计算 μ_i 与 dμ_i/dn_k
  phaseCtx.mu = thermo_.chemicalPotentials(st);
  auto dmun   = thermo_.dmu_dn(st);  // C×C

  // 3) 计算总摩尔数 β_j，并更新 x
  const double beta = std::accumulate(st.moleNumbers.begin(),
                                      st.moleNumbers.end(), 0.0);

  phaseCtx.x.assign(C, 0.0);
  if (beta > 0.0) {
    const double invBeta = 1.0 / beta;
    for (size_t i = 0; i < C; ++i) {
      phaseCtx.x[i] = st.moleNumbers[i] * invBeta;
    }
  } else {
    const double uniform = 1.0 / static_cast<double>(C);
    std::fill(phaseCtx.x.begin(), phaseCtx.x.end(), uniform);
  }

  // 4) 组装 m[i][k] = β_j * ( (1/RT) * dμ_i/dn_k + 1 )
  phaseCtx.m.assign(C, std::vector<double>(C, 0.0));
  for (size_t i = 0; i < C; ++i) {
    for (size_t k = 0; k < C; ++k) {
      phaseCtx.m[i][k] = beta * (dmun[i][k] / RT ) + 1.0;
    }
  }

  // 5) 把可能被裁剪过的 n 写回 phaseCtx.state
  phaseCtx.state.moleNumbers = st.moleNumbers;
}

void RandFlash::fixPhaseHessian(PhaseContext& phaseCtx,
                                const PhaseFixOptions& opt)
{
  const size_t C = phaseCtx.state.moleNumbers.size();
  if (C == 0) {
    phaseCtx.m.clear();
    phaseCtx.M.clear();
    phaseCtx.x.clear();
    phaseCtx.mu.clear();
    return;
  }

  // 正常情况下，updatePhaseChemistry 已经算好了 x；
  // 万一没算，兜底成均匀分布
  if (phaseCtx.x.size() != C) {
    phaseCtx.x.assign(C, 1.0 / static_cast<double>(C));
  }

  // 1) 相内 Hessian SPD 修正
  PhaseFixResult fixRes =
      fix_phase_hessian_one_phase(linearSolver_, phaseCtx.x, phaseCtx.m, opt);

  if (fixRes.applied) {
    std::cout << "[PhaseFix] lam_min "
              << fixRes.lam_min_before << " -> "
              << fixRes.lam_min_after << std::endl;
  }

  // 2) 用返回的 m_fixed（如果有），否则就用当前 m
  if (!fixRes.m_fixed.empty()) {
    phaseCtx.m = fixRes.m_fixed;
  }
  // 理论上 fixRes.m_fixed 应该总是非空；加个兜底也无妨
  if (phaseCtx.m.empty()) {
    phaseCtx.m.assign(C, std::vector<double>(C, 0.0));
  }

  // 3) 每次迭代都基于最新的 m 重新计算 M
  if (!fixRes.M_fixed.empty()) {
    // 如果 PhaseFix 内部已经顺便求了逆，就直接用它
    phaseCtx.M = fixRes.M_fixed;
  } else {
    // 否则就在这里求一次 SPD 逆
    phaseCtx.M = ls::invert(phaseCtx.m, linearSolver_);
  }
}

void randflash::RandFlash::assembleLocalJacobian(
  const thermo::PhaseState& state,
  std::vector<std::vector<double>>& m,
  std::vector<double>& mu)
{
  PhaseContext ctx;
  ctx.state = state;

  updatePhaseChemistry(ctx);

  m  = ctx.m;
  mu = ctx.mu;
}


// 2) 全局系统装配：对应论文式 (4.21)-(4.25)
//    构造 (E+2)x(E+2) 系数矩阵 Acoef 和 RHS rhs
void RandFlash::assembleGlobalSystem(
  double temperature,
  const std::vector<std::vector<std::vector<double>>>& Ms,
  const std::vector<std::vector<double>>& mus,
  const std::vector<std::vector<double>>& nPhases,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& feedComposition, 
  std::vector<double>& Acoef,
  std::vector<double>& rhs)
{
  const int E = static_cast<int>(elementMatrix.size());
  const int C = static_cast<int>(elementMatrix[0].size());
  const int F = static_cast<int>(Ms.size());

  const int N = E + F;
  Acoef.assign(N * N, 0.0);
  rhs.assign(N, 0.0);

  const double RT = R_CONST * temperature;

  // 1) 每个相的 β_j 和 x_j
  std::vector<double> beta(F, 0.0);
  std::vector<std::vector<double>> x(F, std::vector<double>(C, 0.0));
  for (int j = 0; j < F; ++j) {
    const auto& n = nPhases[j];
    beta[j] = std::accumulate(n.begin(), n.end(), 0.0);
    for (int i = 0; i < C; ++i) {
      x[j][i] = (beta[j] > 0.0) ? (n[i] / beta[j]) : 0.0;
    }
  }

  // 2) 左上块：A (∑_j β_j M_j) A^T
  for (int ell = 0; ell < E; ++ell) {
    for (int k = 0; k < E; ++k) {
      double sumjk = 0.0;
      for (int j = 0; j < F; ++j) {
        const auto& Mj = Ms[j];
        const double betaj = beta[j];
        for (int i = 0; i < C; ++i) {
          const double A_li = elementMatrix[ell][i];
          if (A_li == 0.0) continue;
          for (int p = 0; p < C; ++p) {
            const double A_kp = elementMatrix[k][p];
            if (A_kp == 0.0) continue;
            sumjk += A_li * (betaj * Mj[i][p]) * A_kp;
          }
        }
      }
      Acoef[ell * N + k] = sumjk;
    }
  }

  // 3) 左下 / 右上块：A X（每个相一列）
  for (int ell = 0; ell < E; ++ell) {
    for (int j = 0; j < F; ++j) {
      double AX = 0.0;
      for (int i = 0; i < C; ++i) {
        AX += elementMatrix[ell][i] * x[j][i];
      }
      const int col = E + j;
      Acoef[ell * N + col] = AX;
      Acoef[col * N + ell] = AX;  // 对称
    }
  }

  // ===== RHS 部分 =====

  // 4) u1: 
  for (int ell = 0; ell < E; ++ell) {
    double val = 0.0;
    
    // Term 1: Potential contribution
    for (int j = 0; j < F; ++j) {
      const auto& Mj  = Ms[j];
      const auto& muj = mus[j];
      const double betaj = beta[j];
      for (int i = 0; i < C; ++i) {
        const double A_li = elementMatrix[ell][i];
        if (A_li == 0.0) continue;
        for (int p = 0; p < C; ++p) {
          val += A_li * (betaj * Mj[i][p] * (muj[p] / RT));
        }
      }
    }

    rhs[ell] = val;
  }

  // 5) u2: 每个相的“还原自由能”部分 (G_j/RT)
  for (int j = 0; j < F; ++j) {
    const auto& n   = nPhases[j];
    const auto& muj = mus[j];
    const double betaj = beta[j];
    double red = 0.0;
    if (betaj > 0.0) {
      for (int i = 0; i < C; ++i) {
        red += (n[i] / betaj) * (muj[i] / RT);
      }
    }
    rhs[E + j] = red;
  }
}

// 3) 线搜索：对应论文式 (4.32)
double RandFlash::lineSearch(
  const std::vector<std::vector<double>>& nPhases,
  const std::vector<std::vector<double>>& dnPhases,
  const std::vector<std::vector<double>>& gPhases, // mu/RT
  double init_alpha,
  double min_alpha,
  double shrink
) {
  auto dot = [](const std::vector<double>& a, const std::vector<double>& b){
      double s=0.0; for (size_t i=0;i<a.size();++i) s += a[i]*b[i]; return s;
  };
  
  const size_t F = nPhases.size();
  if (F == 0) return min_alpha;
  const size_t C = nPhases[0].size();

  // 检查 alpha 步长后是否保持正值
  const double significant_mole = 1e-10;

  auto positive_after = [&](double a){
      for (size_t j=0; j<F; ++j) {
          for (size_t i=0; i<C; ++i) {
              if (nPhases[j][i] > significant_mole) {
                if (nPhases[j][i] + a*dnPhases[j][i] <= 0.0) return false;
              }
          }
      }
      return true;
  };

  // 下降方向度量: sum_j (dn_j . g_j)
  double dir = 0.0;
  for (size_t j=0; j<F; ++j) {
      dir += dot(dnPhases[j], gPhases[j]);
  }

  const double dir_eps = 1e-12; 
  std::cout << "  descent metric = " << dir << std::endl; 

  double alpha = init_alpha;
  while (alpha > min_alpha) {
    bool positive = positive_after(alpha);
    bool descent  = (dir < -dir_eps) ? (alpha * dir < 0.0)
                                     : (std::abs(dir) <= dir_eps); 
    if (positive) break;
    alpha *= shrink;
  }
  if (alpha <= min_alpha) alpha = min_alpha;
  return alpha;
}

// alpha下限可变版的线搜索：先做可行性裁剪，再检查下降性 （暂无使用）
static double lineSearchFeasible(
  const std::vector<double>& nV, const std::vector<double>& nL,
  const std::vector<double>& dnV, const std::vector<double>& dnL,
  const std::vector<double>& gV,  const std::vector<double>& gL,
  double init_alpha, double min_alpha, double shrink)
{
  auto maxFeasibleAlpha1 = [](const std::vector<double>& n,
                              const std::vector<double>& dn){
      double a = 1.0;
      for (size_t i = 0; i < n.size(); ++i) {
          if (dn[i] < 0.0) {
              // 0.99 给一点余量，避免数值触边
              a = std::min(a, 0.99 * n[i] / (-dn[i]));
          }
      }
      return a;
  };

  // 方向是否下降（对约束自由能）：与 alpha 无关，判一次即可
  auto dirDot = [&](double scale)->double {
      (void)scale; // 方向不依赖 alpha，保持接口一致
      double s = 0.0;
      for (size_t i=0;i<gV.size();++i) s += dnV[i]*gV[i];
      for (size_t i=0;i<gL.size();++i) s += dnL[i]*gL[i];
      return s;
  };
  const double dir = dirDot(1.0);
  const double dir_eps = 1e-12;
  std::cout << "  descent metric = " << dir << std::endl; 
  double alpha = init_alpha; // 仍然以 1.0 起步最稳
  while (alpha > min_alpha) {
      // 先做“可行性裁剪”，再检查下降性
      double aV = maxFeasibleAlpha1(nV, dnV);
      double aL = maxFeasibleAlpha1(nL, dnL);
      alpha = std::min(alpha, std::min(aV, aL));

      if (alpha <= min_alpha) break;

      bool descent = (dir < -dir_eps) ? true : (std::abs(dir) <= dir_eps);
      if (descent) return alpha;

      alpha *= shrink; // 只在不下降时缩步
  }
  return min_alpha;
}

// 4) 解线性系统：当前的 SVD 解法与残差提示
std::vector<double> RandFlash::solveGlobalLinearSystem(
  const std::vector<double>& Acoef,
  const std::vector<double>& rhs,
  double* residual_out) const
{
  const int N = static_cast<int>(rhs.size());
  assert(static_cast<int>(Acoef.size()) == N * N);

  double res_norm = 0.0;
  std::vector<double> sol = linearSolver_.solveDense(N, Acoef, rhs, &res_norm);

  if (residual_out) {
    *residual_out = res_norm;
  }
  if (res_norm > 1e-6) {
    std::cerr << "[RandFlash] Warning: high residual = " << res_norm << "\n";
  }
  return sol;
}

// 5) 结果回代：由 {Λ, Δβ} 得到 Δn；
void RandFlash::backSubstituteDeltas(
  double temperature,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<std::vector<std::vector<double>>>& Ms,
  const std::vector<std::vector<double>>& mus,
  const std::vector<std::vector<double>>& nPhases,
  const std::vector<double>& Lambda,
  const std::vector<double>& deltaBeta,
  std::vector<std::vector<double>>& dnPhases) const
{
  const size_t F = Ms.size();
  const size_t C = nPhases.empty() ? 0 : nPhases[0].size();
  const size_t E = elementMatrix.size();

  dnPhases.assign(F, std::vector<double>(C, 0.0));

  const double RT = R_CONST * temperature;

  // 1) A^T * Lambda, 维度 C
  std::vector<double> lambdaComponent(C, 0.0);
  for (size_t p = 0; p < C; ++p) {
    for (size_t ell = 0; ell < E; ++ell) {
      lambdaComponent[p] += elementMatrix[ell][p] * Lambda[ell];
    }
  }

  // 2) 针对每个相 j 计算 diff_j = A^T λ − μ_j/RT，并回代 Δn^{(j)}
  for (size_t j = 0; j < F; ++j) {
    const auto& Mj  = Ms[j];
    const auto& muj = mus[j];
    const auto& n   = nPhases[j];

    double betaj = std::accumulate(n.begin(), n.end(), 0.0);
    if (betaj <= 0.0) continue;

    std::vector<double> diff(C, 0.0);
    for (size_t p = 0; p < C; ++p) {
      diff[p] = lambdaComponent[p] - muj[p] / RT;
    }

    auto& dn = dnPhases[j];

    for (size_t i = 0; i < C; ++i) {
      const double xij = n[i] / betaj;
      double comb = 0.0;
      for (size_t p = 0; p < C; ++p) {
        comb += Mj[i][p] * diff[p];
      }
      dn[i] = xij * deltaBeta[j] + betaj * comb;
    }
  }

  // 通用打印逻辑，支持任意相数
  // std::cout << "\n  [BackSub] Deltas (dn):";
  // for (size_t j = 0; j < F; ++j) {
  //     std::cout << "\n    Phase " << j << ": ";
  //     for (double val : dnPhases[j]) {
  //         std::cout << val << " ";
  //     }
  // }
  // std::cout << "\n" << std::endl;
}


// 6) 结果更新：线搜索 + n 的更新与“更新后”打印（逻辑不变）
void RandFlash::applyUpdate(
  double temperature,
  const std::vector<std::vector<double>>& mus,
  const std::vector<std::vector<double>>& dnPhases,
  std::vector<std::vector<double>>& nPhases_inout) 
{
  const double RT = R_CONST * temperature;
  const size_t F = nPhases_inout.size();

  // 构造 g = mu/RT
  std::vector<std::vector<double>> gPhases(F);
  for(size_t j=0; j<F; ++j) {
      gPhases[j].resize(mus[j].size());
      for(size_t i=0; i<mus[j].size(); ++i) {
          gPhases[j][i] = mus[j][i] / RT;
      }
  }

  double alpha = lineSearch(
      nPhases_inout, dnPhases, gPhases,
      1.0, 1e-10, 0.5);

  std::cout << "  线搜索步长 alpha = " << std::fixed << std::setprecision(4)<<alpha << std::endl;

  // 更新 n，并执行 Clamping
  const double min_mole_floor = 1e-20; // 物理下限
  for(size_t j=0; j<nPhases_inout.size(); ++j) {
      for(size_t i=0; i<nPhases_inout[j].size(); ++i) {
          double new_val = nPhases_inout[j][i] + alpha * dnPhases[j][i];
          
          // 如果更新后变为负数或极小值，强制拉回下限
          if (new_val < min_mole_floor) {
              new_val = min_mole_floor;
          }
          nPhases_inout[j][i] = new_val;
      }
  }

  // 打印各相总摩尔数 (Beta)
  std::cout << "  更新后各相总摩尔数 (Beta): ";
  for(size_t j=0; j<F; ++j) {
      double beta_j = std::accumulate(nPhases_inout[j].begin(), nPhases_inout[j].end(), 0.0);
      std::cout << std::fixed << std::setprecision(6) << beta_j << " ";
  }
  std::cout << "\n";

  // ========== 新增：打印各相摩尔分率（固定小数位，无科学计数法） ==========
  std::cout << "  更新后各相摩尔分率:\n";
  for(size_t j=0; j<F; ++j) {
      // 计算当前相的总摩尔数
      double beta_j = std::accumulate(nPhases_inout[j].begin(), nPhases_inout[j].end(), 0.0);
      // 防止除以0（理论上不会触发，因有min_mole_floor）
      beta_j = std::max(beta_j, min_mole_floor);

      std::cout << "    相 " << j+1 << ": ";
      for(size_t i=0; i<nPhases_inout[j].size(); ++i) {
          // 计算摩尔分率
          double mole_frac = nPhases_inout[j][i] / beta_j;
          // 输出格式：固定8位小数（无科学计数法），覆盖0~1范围的精度需求
          std::cout << std::fixed << std::setprecision(8) << mole_frac << " ";
      }
      std::cout << "\n";
  }
  std::cout << "\n";
}

// 7) 收敛判断：打印与返回布尔值（阈值保持原样：elem_error < 1e-8）
RandFlash::ConvergenceInfo RandFlash::checkConvergence(
  int iternumber,
  double tol,
  double temperature,
  const std::vector<std::vector<double>>& mus,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<std::vector<double>>& nPhases,
  const std::vector<double>& feedComposition) const
{
  RandFlash::ConvergenceInfo info{};
  const double RT = R_CONST * temperature;
  const size_t F = mus.size();
  const size_t C = feedComposition.size();
  const size_t E = elementMatrix.size();

  // 1) 化学势平衡判据: 计算所有相两两之间的最大无量纲化学势差
  //    或者简单点：计算所有相与第0相的差 (假设第0相存在且稳定)
  //    更稳健的方法：找出当前存在的相（beta > 0），计算它们之间的最大差
  double max_mu_diff = 0.0;
  
  if (F > 1) {
    for (size_t j = 1; j < F; ++j) {
        for (size_t i = 0; i < C; ++i) {
              double diff = std::fabs(mus[j][i] - mus[0][i]) / RT;
              max_mu_diff = std::max(max_mu_diff, diff);
        }
    }
}

  // 2) 元素守恒判据
  // sum_phases (n) - feed
  double elem_err = 0.0;
  for (size_t ell = 0; ell < E; ++ell) {
      double sum_n_ell = 0.0;
      for (size_t j=0; j<F; ++j) {
          // nPhases[j] dot elementMatrix[ell]
          for (size_t i=0; i<C; ++i) {
              sum_n_ell += elementMatrix[ell][i] * nPhases[j][i];
          }
      }
      
      // feed dot elementMatrix[ell] (如果 A 是单位阵，这步简化为 sum_n - feed)
      // 注意：elementMatrix * (sum_n - feed)
      double sum_feed_ell = 0.0;
      for (size_t i=0; i<C; ++i) {
           sum_feed_ell += elementMatrix[ell][i] * feedComposition[i];
      }

      elem_err = std::max(elem_err, std::abs(sum_n_ell - sum_feed_ell));
  }

  std::cout << "Iter "<< iternumber <<" | max_mu_diff (dimless): " << max_mu_diff 
            << " | elem_err: " << elem_err << std::endl;

  info.max_mu_diff = max_mu_diff;
  info.elem_error  = elem_err;
  info.converged   = (max_mu_diff < tol && elem_err < 1e-6);
  return info;
}

void RandFlash::printResult(const MultiFlashResult& res) const {
  // 1. 自动从 Backend 获取组分名称
  std::vector<std::string> compNames = thermo_.getComponentNames();

  std::cout << "\n============================== Flash Results ==============================\n";
  std::cout << std::fixed << std::setprecision(2);
  std::cout << "Temperature: " << res.temperature << " K\n";
  std::cout << "Pressure:    " << res.pressure << " Pa\n";
  std::cout << "Iterations:  " << res.iterations << "\n";
  std::cout << "Status:      " << (res.success ? "CONVERGED" : "FAILED") << "\n";
  std::cout << "Error Norm:  " << std::scientific << std::setprecision(4) << res.mu_infinity_norm << "\n";

  size_t F = res.beta.size();
  if (F == 0) {
      std::cout << "No phases returned.\n";
      std::cout << "===========================================================================\n";
      return;
  }

  size_t C = 0;
  if (!res.n_phase.empty()) C = res.n_phase[0].size();

  // 校验一下组分数量是否匹配，防止越界打印
  if (compNames.size() != C) {
      // 如果数量不对（极少情况），补全或截断，防止 crash
      compNames.resize(C, "Unknown");
  }

  double total_moles = 0.0;
  for (double b : res.beta) total_moles += b;

  std::cout << std::fixed << std::setprecision(5);

  for (size_t j = 0; j < F; ++j) {
      double phase_frac = (total_moles > 1e-12) ? (res.beta[j] / total_moles) : 0.0;
      
      std::string phaseType = "Phase";
      if (j == 0) phaseType = "Vapor (Approx)";
      else phaseType = "Liquid " + std::to_string(j);

      std::cout << "\n---------------------------------------------------------------------------\n";
      std::cout << " " << phaseType << " " << j << " | Phase Fraction (Beta): " << phase_frac << " | Total Moles: " << res.beta[j] << "\n";
      std::cout << "---------------------------------------------------------------------------\n";
      std::cout << "  Idx | " << std::left << std::setw(15) << "Component" << " | " 
                << std::right << std::setw(12) << "Mole Frac (x)" << " | " 
                << std::setw(12) << "Moles (n)" << "\n";
      std::cout << "------+-----------------+--------------+--------------\n";

      for (size_t i = 0; i < C; ++i) {
          double n_i = res.n_phase[j][i];
          double x_i = (res.beta[j] > 1e-20) ? (n_i / res.beta[j]) : 0.0;
          
          std::cout << "  " << std::setw(3) << i << " | " 
                    << std::left << std::setw(15) << compNames[i] << " | " 
                    << std::right << std::setw(12) << x_i << " | " 
                    << std::scientific << std::setprecision(4) << n_i << std::fixed << std::setprecision(5) << "\n";
      }
  }
  std::cout << "===========================================================================\n\n";
}

// randflash::FlashResult RandFlash::solveTwoPhase(
//   const FlashInput& input,
//   const std::vector<std::vector<double>>& elementMatrix,
//   const std::vector<double>& initialVaporComposition,
//   const std::vector<double>& initialLiquidComposition,
//   int maxIter,
//   double tol)
// {
//     // 1. 内部构建 SystemContext (对用户隐藏)
//     SystemContext sys;
//     sys.temperature = input.temperature;
//     sys.pressure    = input.pressure;
//     sys.feedMoles   = input.feedMoles;
//     sys.elementMatrix = elementMatrix;
//     sys.phases.resize(2); // 强制两相

//     // 2. 调用通用的初始化逻辑
//     auto init = initializeTwoPhase(sys, initialVaporComposition, initialLiquidComposition);
    
//     // 3. 将初始化结果填入 SystemContext
//     //    约定: Phase 0 = Vapor, Phase 1 = Liquid
//     sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], thermo_.vaporPhaseFlag()};
//     sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], thermo_.liquidPhaseFlag()};

//     // 4. 委托给核心求解器
//     MultiFlashResult multiRes = solveGeneral(sys, maxIter, tol);

//     // 5. 结果适配回旧的 FlashResult
//     FlashResult res;
//     res.success = multiRes.success;
//     res.pressure = multiRes.pressure;
//     res.temperature = multiRes.temperature;
//     res.feedComposition = multiRes.feedComposition;
//     res.iterations = multiRes.iterations;
//     res.convergenceError = multiRes.mu_infinity_norm;

//     if (multiRes.success && multiRes.n_phase.size() >= 2) {
//         // 映射回 vapor/liquid
//         res.vaporComposition = multiRes.n_phase[0];
//         res.liquidComposition = multiRes.n_phase[1];
        
//         double bV = multiRes.beta[0];
//         double bL = multiRes.beta[1];
//         double total = bV + bL;
//         res.vaporFraction = (total > 0.0) ? (bV / total) : 0.0;
//     }

//     return res;
// }

