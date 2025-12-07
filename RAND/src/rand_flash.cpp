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


// === 在切空间做 SPD 修正，并保持 m x = 1 ===
static PhaseFixResult fix_phase_hessian_one_phase(
  ls::LinearSolverInterface& solver,
  const std::vector<double>& x,                 // 相内 x
  const std::vector<std::vector<double>>& m_in, // 相内 m
  const PhaseFixOptions& opt = {})
{
  const int C = static_cast<int>(x.size());
  if (C == 0) {
    PhaseFixResult out;
    out.applied = false;
    out.lam_min_before = 0.0;
    out.lam_min_after  = 0.0;
    out.m_fixed = m_in;
    out.M_fixed = m_in;
    return out;
  }
  if (static_cast<int>(m_in.size()) != C) {
    throw std::invalid_argument("fix_phase_hessian_one_phase: m_in row size mismatch");
  }
  for (int i = 0; i < C; ++i) {
    if (static_cast<int>(m_in[i].size()) != C) {
      throw std::invalid_argument("fix_phase_hessian_one_phase: m_in must be CxC");
    }
  }

  // 1) 对称部 ms = 0.5 * (m + m^T)
  std::vector<std::vector<double>> ms(C, std::vector<double>(C, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int j = 0; j < C; ++j) {
      ms[i][j] = 0.5 * (m_in[i][j] + m_in[j][i]);
    }
  }

  // 2) 切空间基 B 和切空间内 Hessian mt = B^T ms B
  auto B = ls::tangentBasis(C);     // C × (C-1)
  const int T = C - 1;

  // tmp = ms * B   (C × T)
  std::vector<std::vector<double>> tmp(C, std::vector<double>(T, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int k = 0; k < T; ++k) {
      double s = 0.0;
      for (int j = 0; j < C; ++j) {
        s += ms[i][j] * B[j][k];
      }
      tmp[i][k] = s;
    }
  }

  // mt = B^T * tmp  (T × T)
  std::vector<std::vector<double>> mt(T, std::vector<double>(T, 0.0));
  for (int p = 0; p < T; ++p) {
    for (int q = 0; q < T; ++q) {
      double s = 0.0;
      for (int i = 0; i < C; ++i) {
        s += B[i][p] * tmp[i][q];
      }
      mt[p][q] = s;
    }
  }

  // 3) 计算切空间最小特征值
  std::vector<double> mt_flat(T * T);
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < T; ++j) {
      mt_flat[i * T + j] = mt[i][j];
    }
  }
  std::vector<double> evals;
  std::vector<double> evecs;
  solver.eigenDecomposeSymmetric(T, mt_flat, evals, evecs);

  PhaseFixResult out;
  if (evals.empty()) {
    out.applied = false;
    out.lam_min_before = 0.0;
    out.lam_min_after  = 0.0;
    out.m_fixed = m_in;
    return out;
  }

  const double lam_min = evals[0]; // 升序
  out.lam_min_before = lam_min;
  // std::cout<<"lam_min: "<<lam_min;

  double add = 0.0;
  std::vector<std::vector<double>> mt_fixed = mt;
  if (lam_min < opt.eig_floor) {
    add = opt.eig_floor - lam_min + opt.eps_shift;
    for (int i = 0; i < T; ++i) {
      mt_fixed[i][i] += add;
    }
  }

  if (lam_min >= opt.eig_floor) {
    // 不需要修正
    out.applied = false;
    out.lam_min_after = lam_min;
    out.m_fixed = m_in;
    return out;
  }

  // 4) 只回写切空间修正：m_tan = B * mt_fixed * B^T
  std::vector<std::vector<double>> tmp2(C, std::vector<double>(T, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int k = 0; k < T; ++k) {
      double s = 0.0;
      for (int p = 0; p < T; ++p) {
        s += B[i][p] * mt_fixed[p][k];
      }
      tmp2[i][k] = s;
    }
  }

  std::vector<std::vector<double>> m_tan(C, std::vector<double>(C, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int j = 0; j < C; ++j) {
      double s = 0.0;
      for (int k = 0; k < T; ++k) {
        s += tmp2[i][k] * B[j][k];
      }
      m_tan[i][j] = s;
    }
  }

  // 5) 构造 K = 1 z^T + z 1^T + gamma 11^T，强制 m x = 1
  std::vector<double> ones(C, 1.0);

  // r = 1 - m_tan * x
  std::vector<double> mx(C, 0.0);
  for (int i = 0; i < C; ++i) {
    double s = 0.0;
    for (int j = 0; j < C; ++j) {
      s += m_tan[i][j] * x[j];
    }
    mx[i] = s;
  }

  std::vector<double> r(C);
  for (int i = 0; i < C; ++i) {
    r[i] = 1.0 - mx[i];
  }

  // 分解 r = r_perp + rho * 1
  double rho = 0.0;
  for (int i = 0; i < C; ++i) rho += r[i];
  rho /= static_cast<double>(C);

  std::vector<double> r_perp(C);
  for (int i = 0; i < C; ++i) {
    r_perp[i] = r[i] - rho * ones[i]; // ones[i]==1
  }

  // 设 z = r_perp, gamma = rho - r_perp^T x
  std::vector<double> z = r_perp;
  double rperp_dot_x = 0.0;
  for (int i = 0; i < C; ++i) {
    rperp_dot_x += r_perp[i] * x[i];
  }
  double gamma = rho - rperp_dot_x;

  // K = 1 z^T + z 1^T + gamma 11^T
  std::vector<std::vector<double>> m_fix(C, std::vector<double>(C, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int j = 0; j < C; ++j) {
      double Kij = ones[i] * z[j] + z[i] * ones[j] + gamma * ones[i] * ones[j];
      m_fix[i][j] = m_tan[i][j] + Kij;
    }
  }

  // （可选）数值余量：轻微对称化，避免舍入
  for (int i = 0; i < C; ++i) {
    for (int j = i + 1; j < C; ++j) {
      double avg = 0.5 * (m_fix[i][j] + m_fix[j][i]);
      m_fix[i][j] = m_fix[j][i] = avg;
    }
  }

  // 6) 检查切空间最小特征值（B^T m_fix B）
  std::vector<std::vector<double>> tmp3(C, std::vector<double>(T, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int k = 0; k < T; ++k) {
      double s = 0.0;
      for (int j = 0; j < C; ++j) {
        s += m_fix[i][j] * B[j][k];
      }
      tmp3[i][k] = s;
    }
  }

  std::vector<std::vector<double>> mt_chk(T, std::vector<double>(T, 0.0));
  for (int p = 0; p < T; ++p) {
    for (int q = 0; q < T; ++q) {
      double s = 0.0;
      for (int i = 0; i < C; ++i) {
        s += B[i][p] * tmp3[i][q];
      }
      mt_chk[p][q] = s;
    }
  }

  std::vector<double> mt_chk_flat(T * T);
  for (int i = 0; i < T; ++i) {
    for (int j = 0; j < T; ++j) {
      mt_chk_flat[i * T + j] = mt_chk[i][j];
    }
  }
  std::vector<double> evals_chk;
  std::vector<double> evecs_chk;
  solver.eigenDecomposeSymmetric(T, mt_chk_flat, evals_chk, evecs_chk);
  if (!evals_chk.empty()) {
    out.lam_min_after = evals_chk[0];
  } else {
    out.lam_min_after = 0.0;
  }

  // 7) 求逆：M_fix = m_fix^{-1}
  std::vector<double> m_fix_flat(C * C);
  for (int i = 0; i < C; ++i) {
    for (int j = 0; j < C; ++j) {
      m_fix_flat[i * C + j] = m_fix[i][j];
    }
  }
  std::vector<double> M_flat = solver.invertSPD(C, m_fix_flat);
  std::vector<std::vector<double>> M_fix(C, std::vector<double>(C, 0.0));
  for (int i = 0; i < C; ++i) {
    for (int j = 0; j < C; ++j) {
      M_fix[i][j] = M_flat[i * C + j];
    }
  }

  // 不变量核验：M*1 是否等于 x（理论上应当精确成立）
  std::vector<double> M1(C, 0.0);
  for (int i = 0; i < C; ++i) {
    double s = 0.0;
    for (int j = 0; j < C; ++j) {
      s += M_fix[i][j] * ones[j];
    }
    M1[i] = s;
  }

  std::vector<double> diff(C, 0.0);
  double diff2 = 0.0;
  for (int i = 0; i < C; ++i) {
    diff[i] = M1[i] - x[i];
    diff2 += diff[i] * diff[i];
  }

  if (std::sqrt(diff2) > 1e-9) {
    const double invC = 1.0 / static_cast<double>(C);
    // 极小对称化补偿（理论上用不到）
    for (int i = 0; i < C; ++i) {
      for (int j = 0; j < C; ++j) {
        M_fix[i][j] -= 0.5 * (diff[i] + diff[j]) * invC;
      }
    }
  }

  // 输出
  out.applied = (add > 0.0);     // 只有切空间真做了移位才算“应用”
  out.m_fixed = std::move(m_fix);
  out.M_fixed = std::move(M_fix);
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
  std::cout << "\n  [BackSub] Deltas (dn):";
  for (size_t j = 0; j < F; ++j) {
      std::cout << "\n    Phase " << j << ": ";
      for (double val : dnPhases[j]) {
          std::cout << val << " ";
      }
  }
  std::cout << "\n" << std::endl;
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

  std::cout << "  线搜索步长 alpha = " << alpha << std::endl;

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

  // 打印简单的更新信息 (可选: 仅打印 beta 以简化日志)
  std::cout << "  更新后 Betas: ";
  for(size_t j=0; j<F; ++j) {
      double b = std::accumulate(nPhases_inout[j].begin(), nPhases_inout[j].end(), 0.0);
      std::cout << b << " ";
  }
  std::cout << "\n\n";
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

