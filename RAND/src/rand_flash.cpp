#include "rand_flash.hpp"
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
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

  // 工具：归一化/基本检查
  // 替换 initializeTwoPhase 里原来的 normalized(...)
  auto normalized = [](std::vector<double> v){
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    if (!(s > 0.0)) throw std::invalid_argument("Composition sum must be positive");
    for (double xi : v)
      if (xi < -1e-12) throw std::invalid_argument("Composition has negative entry");
    // 轻度纠偏：把极小的负数夹到 0，再归一
    for (double &xi : v) if (xi < 0.0) xi = 0.0;
    double s2 = std::accumulate(v.begin(), v.end(), 0.0);
    if (std::abs(s2 - 1.0) > 1e-8) {
      std::cerr << "[init] composition not normalized, sum=" << s2 << "\n";
    }
    for (double &x : v) x /= (s2 > 0.0 ? s2 : 1.0);
    return v;
  };

  auto near_zero = [](double x){ return std::abs(x) < 1e-14; };

  // 进料分率 z
  std::vector<double> z(C);
  for (size_t i=0;i<C;++i) z[i] = sys.feedMoles[i] / Ftot;

  // 结果
  InitResult out;

  const double eps = 1e-14;
  const double Fmin = eps, Fmax = Ftot * (1.0 - margin);

  // CASE A：只给了气相组成 y
  if (!vaporGuess.empty() && liquidGuess.empty()) {
    if (vaporGuess.size() != C) throw std::invalid_argument("vaporGuess size mismatch");
    out.y0 = normalized(vaporGuess);

    // 可行域上界：保证 x_i = (F_i - beta*y_i)/(Ftot - beta) >= 0
    double ub = Fmax; // 也要避免 beta→Ftot
    for (size_t i=0;i<C;++i) {
      if (out.y0[i] > eps) {
        ub = std::min(ub, sys.feedMoles[i] / out.y0[i]);
      }
    }
    ub = std::max(ub, Fmin);             // 上界至少要大于 0
    double beta = std::min(0.5*Ftot, ub*(1.0 - margin));  // 取居中且离边界有裕度

    // 由 y 和 beta 反解 x，并保证分母安全
    const double Ltot = Ftot - beta;
    if (Ltot <= Ftot*margin) beta = Ftot*(1.0 - margin);
    std::vector<double> x(C);
    for (size_t i=0;i<C;++i) {
      x[i] = (sys.feedMoles[i] - beta*out.y0[i]) / (Ftot - beta);
      if (x[i] < 0.0) x[i] = 0.0; //（理论上不会触发，容错）
    }
    // 轻度重归一，消除舍入误差
    x = normalized(x);
    out.x0 = x;

    out.betaV = beta;
    out.betaL = Ftot - beta;
    out.nV.resize(C); out.nL.resize(C);
    for (size_t i=0;i<C;++i) {
      out.nV[i] = out.betaV * out.y0[i];
      out.nL[i] = out.betaL * out.x0[i];
    }
    return out;
  }

  // CASE B：只给了液相组成 x
  if (vaporGuess.empty() && !liquidGuess.empty()) {
    if (liquidGuess.size() != C) throw std::invalid_argument("liquidGuess size mismatch");
    out.x0 = normalized(liquidGuess);

    // 可行域下界：保证 y_i = (F_i - (Ftot - beta)*x_i)/beta >= 0
    // 即 F_i - (Ftot - beta)*x_i >= 0  =>  beta >= Ftot - F_i/x_i  (x_i>0 时)
    double lb = Fmin;
    for (size_t i=0;i<C;++i) {
      if (out.x0[i] > eps) {
        lb = std::max(lb, Ftot - sys.feedMoles[i]/out.x0[i]);
      }
    }
    lb = std::max(lb, Fmin);
    double beta = std::max(0.5*Ftot, lb*(1.0 + margin));
    if (beta >= Fmax) beta = 0.5*(lb + Fmax); // 夹到可行区间内

    // 由 x 和 beta 反解 y
    std::vector<double> y(C);
    for (size_t i=0;i<C;++i) {
      y[i] = (sys.feedMoles[i] - (Ftot - beta)*out.x0[i]) / beta;
      if (y[i] < 0.0) y[i] = 0.0; //（理论上不会触发，容错）
    }
    y = normalized(y);
    out.y0 = y;

    out.betaV = beta;
    out.betaL = Ftot - beta;
    out.nV.resize(C); out.nL.resize(C);
    for (size_t i=0;i<C;++i) {
      out.nV[i] = out.betaV * out.y0[i];
      out.nL[i] = out.betaL * out.x0[i];
    }
    return out;
  }

  // CASE C：同时给了 y 与 x —— 用最小二乘估计“共识” beta，再恢复 nV/nL
  if (!vaporGuess.empty() && !liquidGuess.empty()) {
    if (vaporGuess.size() != C || liquidGuess.size() != C)
      throw std::invalid_argument("Initial composition size mismatch");
    out.y0 = normalized(vaporGuess);
    out.x0 = normalized(liquidGuess);

    // 目标：最小化 ∑_i [ beta*(y_i - x_i) - Ftot*(z_i - x_i) ]^2
    double num = 0.0, den = 0.0;
    for (size_t i=0;i<C;++i) {
      const double d = out.y0[i] - out.x0[i];
      num += d * (z[i] - out.x0[i]);
      den += d * d;
    }
    double beta = (den > 0.0) ? Ftot * (num / den) : 0.5*Ftot;

    // 夹到可行区间：保证两相组成非负且分母不接近 0
    // 上界来自 A 案例，下界来自 B 案例
    double ub = Fmax, lb = Fmin;
    for (size_t i=0;i<C;++i) {
      if (out.y0[i] > eps) ub = std::min(ub, sys.feedMoles[i] / out.y0[i]);
      if (out.x0[i] > eps) lb = std::max(lb, Ftot - sys.feedMoles[i]/out.x0[i]);
    }
    beta = std::min(std::max(beta, lb*(1.0 + margin)), ub*(1.0 - margin));

    out.betaV = beta;
    out.betaL = Ftot - beta;

    out.nV.resize(C); out.nL.resize(C);
    for (size_t i=0;i<C;++i) {
      out.nV[i] = out.betaV * out.y0[i];
      out.nL[i] = out.betaL * out.x0[i];
    }
    return out;
  }

  // CASE D：既没有给 y 也没有给 x —— 用 Thermopack 的 Wilson K 做两相初始化
  {
    // 1) feed mole fraction z 已经在函数开头算过：z[i] = sys.feedMoles[i] / Ftot;

    // 2) Thermopack 计算 Wilson K 值
    std::vector<double> K(C, 0.0);
    {
      double T = sys.temperature; // 或改成传进来的 T
      double P = sys.pressure; // 或改成传进来的 P
      thermo_.wilsonK(T, P, K);
    }

    // 3) 用 Rachford-Rice ∑ z_i (K_i-1)/(1+β(K_i-1)) = 0 解出蒸汽分率 β ∈ [0,1]
    auto rr = [&](double beta) {
      double s = 0.0;
      for (size_t i = 0; i < C; ++i) {
        const double d = K[i] - 1.0;
        double denom = 1.0 + beta * d;
        if (denom < 1e-12) denom = 1e-12; // 数值保护
        s += z[i] * d / denom;
      }
      return s;
    };

    double beta_lo = margin;
    double beta_hi = 1.0 - margin;
    double f_lo = rr(beta_lo);
    double f_hi = rr(beta_hi);

    double beta = 0.5 * (beta_lo + beta_hi);

    // 如果区间两端同号，说明 Wilson 在此 T,P 下不太“像两相”，简单退回 β=0.5
    if (f_lo * f_hi < 0.0) {
      for (int it = 0; it < 40; ++it) {
        beta = 0.5 * (beta_lo + beta_hi);
        double f = rr(beta);
        if (std::fabs(f) < 1e-12) break;
        if (f * f_lo > 0.0) {
          beta_lo = beta;
          f_lo = f;
        } else {
          beta_hi = beta;
          f_hi = f;
        }
      }
    } else {
      beta = 0.5; // 单相附近 → 随便给个 0.5 的两相初值
    }

    // 4) 根据 β 和 K 求 x_i / y_i
    out.x0.assign(C, 0.0);
    out.y0.assign(C, 0.0);
    for (size_t i = 0; i < C; ++i) {
      const double d = K[i] - 1.0;
      double denom = 1.0 + beta * d;
      if (denom < 1e-12) denom = 1e-12;
      out.x0[i] = z[i] / denom;
      out.y0[i] = K[i] * out.x0[i];
      if (out.x0[i] < 0.0) out.x0[i] = 0.0;
      if (out.y0[i] < 0.0) out.y0[i] = 0.0;
    }

    // 轻度归一，消除数值误差
    out.x0 = normalized(out.x0);
    out.y0 = normalized(out.y0);

    // 5) 转成 nV / nL
    out.betaV = beta * Ftot;
    out.betaL = (1.0 - beta) * Ftot;

    out.nV.resize(C);
    out.nL.resize(C);
    for (size_t i = 0; i < C; ++i) {
      out.nV[i] = out.betaV * out.y0[i];
      out.nL[i] = out.betaL * out.x0[i];
    }

    return out;
  }



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

  // 4) u1: 元素守恒相关项 (文中对应式 (4.22) 那个 μ 部分)
  for (int ell = 0; ell < E; ++ell) {
    double val = 0.0;
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

  // 5) u2: 每个相的“还原自由能”部分
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

  // feedComposition 当前仍未显式用到（和你现有实现一致），先保留接口不删
}

// 3) 线搜索：对应论文式 (4.32)
double RandFlash::lineSearch(
  const std::vector<double>& nV,
  const std::vector<double>& nL,
  const std::vector<double>& dnV,
  const std::vector<double>& dnL,
  const std::vector<double>& gV,   // = muV/RT  或 ln(fhatV)
  const std::vector<double>& gL,   // = muL/RT  或 ln(fhatL)
  double init_alpha,               // = 1.0
  double min_alpha,                // = 1e-10
  double shrink                    // = 0.5
) {
  auto dot = [](const std::vector<double>& a, const std::vector<double>& b){
      double s=0.0; for (size_t i=0;i<a.size();++i) s += a[i]*b[i]; return s;
  };
  auto positive_after = [&](double a){
      for (size_t i=0;i<nV.size();++i)
          if (nV[i] + a*dnV[i] <= 0.0 || nL[i] + a*dnL[i] <= 0.0) return false;
      return true;
  };

  // 方向的“下降性”度量：sum_j sum_i Δn_{i,j} * g_{i,j}
  const double dir = dot(dnV, gV) + dot(dnL, gL);
  const double dir_eps = 1e-12; // 允许非常小的非负/非负数噪声
  std::cout << "  descent metric = " << dir << std::endl; // 应该是负的；若非负，就会触发回溯
  double alpha = init_alpha;
  while (alpha > min_alpha) {
    bool positive = positive_after(alpha);
    bool descent  = (dir < -dir_eps) ? (alpha * dir < 0.0)
                                     : (std::abs(dir) <= dir_eps); // 近似驻点也放行
    if (positive && descent) break;
    alpha *= shrink;
  }
  if (alpha <= min_alpha) alpha = min_alpha; // 兜底
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

// 5) 结果回代：由 {Λ, Δβ} 得到 {ΔnV, ΔnL}；保持打印 dnV/dnL
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

  // 为了不破坏你原来的调试输出风格，如果是两相就按原来的格式打印
  if (F == 2) {
    const auto& dnV = dnPhases[0];
    const auto& dnL = dnPhases[1];
    std::cout << "\n  dnV: ";
    for (double dv : dnV) std::cout << dv << " ";
    std::cout << "\n  dnL: ";
    for (double dv : dnL) std::cout << dv << " ";
    std::cout << "\n" << std::endl;
  }
}


// 6) 结果更新：线搜索 + n 的更新与“更新后”打印（逻辑不变）
double RandFlash::applyUpdate(
  double temperature,
  const std::vector<double>& muV,
  const std::vector<double>& muL,
  const std::vector<double>& dnV,
  const std::vector<double>& dnL,
  std::vector<double>& nV_inout,
  std::vector<double>& nL_inout) 
{
  const double RT = R_CONST * temperature;

  // 与原实现一致：用 mu/RT 做“方向度量”
  std::vector<double> muV_over_RT(muV.size());
  std::vector<double> muL_over_RT(muL.size());
  for (size_t i = 0; i < muV.size(); ++i) muV_over_RT[i] = muV[i] / RT;
  for (size_t i = 0; i < muL.size(); ++i) muL_over_RT[i] = muL[i] / RT;

  double alpha = lineSearch(
      nV_inout, nL_inout,
      dnV, dnL,
      muV_over_RT, muL_over_RT,
      /*init_alpha=*/1.0, /*min_alpha=*/1e-10, /*shrink=*/0.5);

  std::cout << "  线搜索步长 alpha = " << alpha << std::endl;

  // 更新 nV/nL
  const size_t C = nV_inout.size();
  for (size_t i = 0; i < C; ++i) {
    nV_inout[i] += alpha * dnV[i];
    nL_inout[i] += alpha * dnL[i];
  }

  // 打印“更新后”与新的 x
  double betaV = std::accumulate(nV_inout.begin(), nV_inout.end(), 0.0);
  double betaL = std::accumulate(nL_inout.begin(), nL_inout.end(), 0.0);
  std::cout << "\n  更新后 : ";
  std::cout << "\n  nV: ";
  for (double v : nV_inout) std::cout << v << " ";
  std::cout << "\n  nL: ";
  for (double v : nL_inout) std::cout << v << " ";
  std::cout << "\n  betaV: " << betaV << ", betaL: " << betaL;
  std::cout << "\n  xV: ";
  for (double v : nV_inout) std::cout << (betaV > 1e-12 ? v / betaV : 0.0) << " ";
  std::cout << "\n  xL: ";
  for (double v : nL_inout) std::cout << (betaL > 1e-12 ? v / betaL : 0.0) << " ";
  std::cout << "\n" << std::endl;

  return alpha;
}

// 7) 收敛判断：打印与返回布尔值（阈值保持原样：elem_error < 1e-8）
RandFlash::ConvergenceInfo RandFlash::checkConvergence(
  int iternumber,
  double tol,
  double temperature,
  const std::vector<double>& muV,
  const std::vector<double>& muL,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& nV,
  const std::vector<double>& nL,
  const std::vector<double>& feedComposition) const
{
  RandFlash::ConvergenceInfo info{};
  int C = static_cast<int>(nV.size());
  int E = static_cast<int>(elementMatrix.size());

  // 1) 无量纲化的 max |(muV - muL)/RT|
  const double RT = R_CONST * temperature;   // 和 assembleLocalJacobian 里保持一致
  double max_mu_diff = 0.0;
  for (size_t i = 0; i < muV.size(); ++i) {
    double dimless_diff = std::fabs(muV[i] - muL[i]) / RT;
    max_mu_diff = std::max(max_mu_diff, dimless_diff);
  }

  // 2) 元素守恒 L∞ 范数（维持原实现）
  auto elementResidualInf = [&](const std::vector<std::vector<double>>& A,
                                const std::vector<double>& nV_,
                                const std::vector<double>& nL_,
                                const std::vector<double>& feed) {
    const size_t E = A.size(), C = feed.size();
    double infnorm = 0.0;
    for (size_t ell = 0; ell < E; ++ell) {
      double acc = 0.0;
      for (size_t i = 0; i < C; ++i) {
        acc += A[ell][i] * ((nV_[i] + nL_[i]) - feed[i]);
      }
      infnorm = std::max(infnorm, std::abs(acc));
    }
    return infnorm;
  };
  double elem_err = elementResidualInf(elementMatrix, nV, nL, feedComposition);

  std::cout << "Iter "<< iternumber <<" | max_mu_diff (dimless): " << max_mu_diff << std::endl;
  std::cout << "Iter "<< iternumber <<" | element error: " << elem_err << "\n" << std::endl;

  info.max_mu_diff = max_mu_diff;   // 现在这个就是无量纲差
  info.elem_error  = elem_err;
  info.converged   = (max_mu_diff < tol && elem_err < 1e-6);
  return info;
}


randflash::FlashResult RandFlash::solveTwoPhase(
  const FlashInput& input,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& initialVaporComposition,   // 可空
  const std::vector<double>& initialLiquidComposition,  // 可空
  int maxIter,
  double tol)
{
  const size_t C = input.feedMoles.size();
  const size_t E = elementMatrix.size();

  // === 0) 构造系统级上下文 ===
  SystemContext sys;
  sys.temperature    = input.temperature;
  sys.pressure       = input.pressure;
  sys.feedMoles      = input.feedMoles;
  sys.elementMatrix  = elementMatrix;
  sys.phases.resize(2);  // 0: vapor, 1: liquid

  // === 1) 初始气液摩尔数 nV, nL ===
  std::vector<double> nV(C), nL(C);
  auto init = initializeTwoPhase(
      sys,
      initialVaporComposition,
      initialLiquidComposition,
      /*margin=*/1e-8);

  nV = init.nV;
  nL = init.nL;

  std::cout << "Calculated betaV: " << init.betaV
            << ", betaL: " << init.betaL << "\n"
            << "初始气相摩尔分率 betaV/(betaV+betaL) = "
            << (init.betaV / (init.betaV + init.betaL)) << std::endl;

  // === 2) 用 PhaseContext 封装两相状态 ===
  auto& vapCtx = sys.phases[0];
  auto& liqCtx = sys.phases[1];

  vapCtx.state = thermo::PhaseState{
    input.temperature, input.pressure, nV, thermo_.vaporPhaseFlag()};
  liqCtx.state = thermo::PhaseState{
    input.temperature, input.pressure, nL, thermo_.liquidPhaseFlag()};

  FlashResult result;
  result.pressure        = input.pressure;
  result.temperature     = input.temperature;
  result.feedComposition = input.feedMoles;
  result.success         = false;

  // 预分配迭代所需容器
  std::vector<double> Acoef, rhs, sol;
  std::vector<double> Lambda;              // Δλ
  std::vector<double> deltaBeta(2);        // Δβ_v, Δβ_l
  std::vector<double> dnV(C), dnL(C);      // Δn^V, Δn^L

  std::cout << "进入迭代循环, maxIter = " << maxIter << std::endl;

  try {
    for (int iter = 0; iter < maxIter; ++iter) {
      std::cout << "\n=== RandFlash iter " << (iter + 1) << " ===\n";

      // --- 2.1 更新 PhaseContext 中的 n ---
      vapCtx.state.moleNumbers = nV;
      liqCtx.state.moleNumbers = nL;

      // --- 2.2 计算局部 μ / m（并更新 x） ---
      updatePhaseChemistry(vapCtx);
      updatePhaseChemistry(liqCtx);

      std::cout << "局部 μ / m 计算完成\n";

      // --- 2.3 相内 Hessian 修正 + 求逆，填充 M ---
      PhaseFixOptions pfx;
      fixPhaseHessian(vapCtx, pfx);
      fixPhaseHessian(liqCtx, pfx);

      std::cout << "相内 Hessian 修正完成\n";

      // --- 3) 组装全局线性系统 A x = rhs ---
      Acoef.clear();
      rhs.clear();

      // 3) 全局系统装配（多相通用接口，目前 F=2）
      std::vector<std::vector<std::vector<double>>> Ms(2);
      Ms[0] = vapCtx.M;
      Ms[1] = liqCtx.M;

      std::vector<std::vector<double>> mus(2);
      mus[0] = vapCtx.mu;
      mus[1] = liqCtx.mu;

      std::vector<std::vector<double>> nPhases(2);
      nPhases[0] = nV;
      nPhases[1] = nL;

      assembleGlobalSystem(
        input.temperature,
        Ms,
        mus,
        nPhases,
        elementMatrix,
        input.feedMoles,   // 进料 n^F（暂时未在公式中使用，但保留接口）
        Acoef,
        rhs
      );

      std::cout << "全局系统装配完成\n";

      // --- 4) 解线性系统，得到 Δλ, Δβ ---
      double lin_resid = 0.0;
      sol = solveGlobalLinearSystem(Acoef, rhs, &lin_resid);

      if (sol.size() != E + 2) {
        throw std::runtime_error("solveGlobalLinearSystem returned wrong size");
      }

      Lambda.assign(sol.begin(), sol.begin() + E);
      deltaBeta[0] = sol[E + 0];  // vapor
      deltaBeta[1] = sol[E + 1];  // liquid

      std::cout << "线性系统求解完成, residual = " << lin_resid << "\n";

      // 5) 回代恢复 Δn（通过多相通用接口）
      std::vector<std::vector<std::vector<double>>> Ms_bs(2);
      Ms_bs[0] = vapCtx.M;
      Ms_bs[1] = liqCtx.M;

      std::vector<std::vector<double>> mus_bs(2);
      mus_bs[0] = vapCtx.mu;
      mus_bs[1] = liqCtx.mu;

      std::vector<std::vector<double>> nPhases_bs(2);
      nPhases_bs[0] = nV;
      nPhases_bs[1] = nL;

      std::vector<std::vector<double>> dnPhases;
      backSubstituteDeltas(
        input.temperature,
        elementMatrix,
        Ms_bs,
        mus_bs,
        nPhases_bs,
        Lambda,
        deltaBeta,
        dnPhases
      );

      // 拆回两相增量，保持后续逻辑不变
      dnV = dnPhases[0];
      dnL = dnPhases[1];


      std::cout << "回代 Δn 完成\n";

      // --- 6) 线搜索 + 更新 nV, nL ---
      double alpha = applyUpdate(
        sys.temperature,
        vapCtx.mu, liqCtx.mu,
        dnV, dnL,
        nV, nL);

      std::cout << "线搜索 alpha = " << alpha << "\n";

      // --- 7) 收敛判断 ---
      auto conv = checkConvergence(
        iter + 1,
        tol,
        sys.temperature,
        vapCtx.mu, liqCtx.mu,
        sys.elementMatrix,
        nV, nL,
        sys.feedMoles);

      if (conv.converged) {
        double betaV_fin = std::accumulate(nV.begin(), nV.end(), 0.0);
        double betaL_fin = std::accumulate(nL.begin(), nL.end(), 0.0);

        result.success           = true;
        result.vaporFraction     = betaV_fin / (betaV_fin + betaL_fin);
        result.vaporComposition  = nV;
        result.liquidComposition = nL;
        result.iterations        = iter + 1;
        result.convergenceError  = conv.max_mu_diff;

        std::cout << "RandFlash 收敛, iter = " << result.iterations
                  << ", vaporFraction = " << result.vaporFraction << "\n";
        return result;
      }
    }
  } catch (const std::exception& e) {
    std::cerr << "[RandFlash] Exception caught during iterations: "
              << e.what() << std::endl;
    result.success = false;
    return result;
  }

  // 若 maxIter 内未收敛
  std::cout << "RandFlash 未在 maxIter 内收敛\n";
  result.success = false;
  return result;
}


