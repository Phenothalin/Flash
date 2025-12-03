#include "rand_flash.hpp"
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
using namespace randflash;
using namespace ls;

// 保证逐分量守恒 + 组分单纯形精确成立 (暂无使用)
inline void enforce_component_balance_and_simplex(
  std::vector<double>& nV, std::vector<double>& nL,
  const std::vector<double>& nFeed, // 全流程不变的进料逐分量摩尔数
  double* betaV_out = nullptr, double* betaL_out = nullptr,
  std::vector<double>* xV_out = nullptr, std::vector<double>* xL_out = nullptr)
{
  const int nc = static_cast<int>(nV.size());
  assert(nc == static_cast<int>(nL.size()) &&
         nc == static_cast<int>(nFeed.size()));

  // 1) 逐分量裁剪，避免数值微负/微超
  for (int i = 0; i < nc; ++i) {
    double nv = nV[i];
    const double nf = nFeed[i];

    // 把 nV 限制在 [0, nF]
    if (nv < 0.0) nv = 0.0;
    if (nv > nf)  nv = nf;
    nV[i] = nv;

    // 由守恒强制 nL
    nL[i] = nf - nV[i];

    // 再保一次非负（理论上不会触发，防卫式）
    if (nL[i] < 0.0) { nV[i] += nL[i]; nL[i] = 0.0; }
    if (nV[i] < 0.0) { nL[i] += nV[i]; nV[i] = 0.0; }
  }

  // 2) 精确计算 beta，并构造单纯形上的 x
  double betaV = 0.0;
  double betaL = 0.0;
  for (int i = 0; i < nc; ++i) {
    betaV += nV[i];
    betaL += nL[i];
  }

  if (betaV_out) *betaV_out = betaV;
  if (betaL_out) *betaL_out = betaL;

  if (xV_out) {
    xV_out->assign(nc, 0.0);
    if (betaV > 0.0) {
      const double invBetaV = 1.0 / betaV;
      for (int i = 0; i < nc; ++i) {
        (*xV_out)[i] = nV[i] * invBetaV;
      }
    }
  }
  if (xL_out) {
    xL_out->assign(nc, 0.0);
    if (betaL > 0.0) {
      const double invBetaL = 1.0 / betaL;
      for (int i = 0; i < nc; ++i) {
        (*xL_out)[i] = nL[i] * invBetaL;
      }
    }
  }
}

// === 构造切空间正交基 B: 1^T y = 0 ===
// 用 (e_i - e_C) 做初基，再 Gram–Schmidt 正交化
static std::vector<std::vector<double>> tangentBasis(int C){
  std::vector<std::vector<double>> B(C, std::vector<double>(C-1, 0.0));

  // 初始基：第 k 列为 e_k - e_C
  for (int k = 0; k < C-1; ++k){ // 列 k
    B[k][k]   =  1.0;
    B[C-1][k] = -1.0;
  }

  // Gram-Schmidt 正交化
  for (int j = 0; j < C-1; ++j){
    // 去除在之前列上的分量
    for (int i = 0; i < j; ++i){
      double proj = 0.0;
      for (int r = 0; r < C; ++r) {
        proj += B[r][i] * B[r][j];
      }
      for (int r = 0; r < C; ++r) {
        B[r][j] -= proj * B[r][i];
      }
    }

    // 归一化
    double nrm2 = 0.0;
    for (int r = 0; r < C; ++r) {
      nrm2 += B[r][j] * B[r][j];
    }
    double nrm = std::sqrt(nrm2);
    if (nrm < 1e-14) { // 退化保护：随机扰动再正交
      for (int r = 0; r < C; ++r) {
        B[r][j] = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
      }
      for (int i = 0; i < j; ++i){
        double proj = 0.0;
        for (int r = 0; r < C; ++r) {
          proj += B[r][i] * B[r][j];
        }
        for (int r = 0; r < C; ++r) {
          B[r][j] -= proj * B[r][i];
        }
      }
      nrm2 = 0.0;
      for (int r = 0; r < C; ++r) {
        nrm2 += B[r][j] * B[r][j];
      }
      nrm = std::sqrt(nrm2);
      if (nrm < 1e-14) {
        throw std::runtime_error("tangentBasis: failed to build non-degenerate basis");
      }
    }
    const double inv_nrm = 1.0 / nrm;
    for (int r = 0; r < C; ++r) {
      B[r][j] *= inv_nrm;
    }
  }

  // 可选：检查 1^T B = 0，这里略
  return B;
}

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
  auto B = tangentBasis(C);     // C × (C-1)
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
  const std::vector<double>& feed,
  const double Temperature,
  const double Pressure,
  const std::vector<double>& vaporGuess,
  const std::vector<double>& liquidGuess,
  double margin) const
{
  const size_t C = feed.size();
  auto sum = [](const std::vector<double>& v){ return std::accumulate(v.begin(), v.end(), 0.0); };
  const double Ftot = sum(feed);
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
  for (size_t i=0;i<C;++i) z[i] = feed[i] / Ftot;

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
        ub = std::min(ub, feed[i] / out.y0[i]);
      }
    }
    ub = std::max(ub, Fmin);             // 上界至少要大于 0
    double beta = std::min(0.5*Ftot, ub*(1.0 - margin));  // 取居中且离边界有裕度

    // 由 y 和 beta 反解 x，并保证分母安全
    const double Ltot = Ftot - beta;
    if (Ltot <= Ftot*margin) beta = Ftot*(1.0 - margin);
    std::vector<double> x(C);
    for (size_t i=0;i<C;++i) {
      x[i] = (feed[i] - beta*out.y0[i]) / (Ftot - beta);
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
        lb = std::max(lb, Ftot - feed[i]/out.x0[i]);
      }
    }
    lb = std::max(lb, Fmin);
    double beta = std::max(0.5*Ftot, lb*(1.0 + margin));
    if (beta >= Fmax) beta = 0.5*(lb + Fmax); // 夹到可行区间内

    // 由 x 和 beta 反解 y
    std::vector<double> y(C);
    for (size_t i=0;i<C;++i) {
      y[i] = (feed[i] - (Ftot - beta)*out.x0[i]) / beta;
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
      if (out.y0[i] > eps) ub = std::min(ub, feed[i] / out.y0[i]);
      if (out.x0[i] > eps) lb = std::max(lb, Ftot - feed[i]/out.x0[i]);
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
    // 1) feed mole fraction z 已经在函数开头算过：z[i] = feed[i] / Ftot;

    // 2) Thermopack 计算 Wilson K 值
    std::vector<double> K(C, 0.0);
    {
      double T = Temperature; // 或改成传进来的 T
      double P = Pressure; // 或改成传进来的 P
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
//    m_j[i][k] = β_j * ( (1/RT) * ∂μ_i/∂n_k + 1 )
//    mu_j 已经由 chemicalPotentials 给出
void randflash::RandFlash::assembleLocalJacobian(
  const thermo::PhaseState& state,
  std::vector<std::vector<double>>& m,
  std::vector<double>& mu)
{
  // 1. 复制一份状态，以便安全地改 n_i
  thermo::PhaseState st = state;
  for (double &ni : st.moleNumbers) {
      if (ni <= 1e-10) ni = 1e-10;
  }

  // 2. 计算 μ 和 ∂μ/∂n
  mu = thermo_.chemicalPotentials(st);
  auto dmun = thermo_.dmu_dn(st);
  // for (size_t i = 0; i < moleNumbers.size(); ++i) {
  //   for (size_t k = 0; k < moleNumbers.size(); ++k) {
  //       std::cout << dmun[i][k] << "  ";  // 同一行打印
  //   }
  //   std::cout << std::endl;  // 换行
  // }    
  // 3. 构造 m_j 矩阵：m_j[i][k] = β_j*(dμ/dn/RT +1),此处 β_j = totalMoles
  double totalMoles = std::accumulate(st.moleNumbers.begin(), st.moleNumbers.end(), 0.0);
  double RT = R_CONST * st.Temperature;
  size_t C = st.moleNumbers.size();
  m.assign(C, std::vector<double>(C,0.0));
  for (size_t i=0;i<C;++i){
    for (size_t k=0;k<C;++k){
      double coeff = dmun[i][k]/RT ;      // (4.11) 中括号内
      m[i][k] = totalMoles * coeff + 1.0;           // 乘以 β_j
    }
  }
  // for(size_t i = 0; i < C; ++i) {
  //   for(size_t k = 0; k < C; ++k) {
  //     std::cout << "m[" << i << "][" << k << "] = "
  //               << m[i][k] << std::endl;
  //   }
  // }
  
  // 轻微对角抖动（相对尺度）
  // double dmean = 0.0; for (size_t i=0;i<C;++i) dmean += std::abs(m[i][i]); dmean = std::max(dmean/C, 1.0);
  // for (size_t i=0;i<C;++i) m[i][i] += 1e-10 * dmean;
}

// 2) 全局系统装配：对应论文式 (4.21)-(4.25)
//    构造 (E+2)x(E+2) 系数矩阵 Acoef 和 RHS rhs
void RandFlash::assembleGlobalSystem(
  double temperature,
  const std::vector<std::vector<double>>& MV,
  const std::vector<std::vector<double>>& ML,
  const std::vector<double>& muV,
  const std::vector<double>& muL,
  const std::vector<double>& moleNumbersV,
  const std::vector<double>& moleNumbersL,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& feedComposition, 
  std::vector<double>& Acoef,
  std::vector<double>& rhs)
{
  int E = (int)elementMatrix.size();
  int C = (int)elementMatrix[0].size();
  int F = 2;
  int N = E + F;
  Acoef.assign(N*N, 0.0);
  rhs.assign(N, 0.0);

  double RT = R_CONST * temperature;
  double betaV = std::accumulate(moleNumbersV.begin(), moleNumbersV.end(), 0.0);
  double betaL = std::accumulate(moleNumbersL.begin(), moleNumbersL.end(), 0.0);
  std::vector<double> xV(C), xL(C);
  for(int i=0;i<C;++i){
    xV[i] = moleNumbersV[i] / betaV;
    xL[i] = moleNumbersL[i] / betaL;
  }

  // 1) 左上块 A (sum_j β_j M_j) A^T
  for(int ell = 0; ell < E; ++ell){
    for(int k = 0; k < E; ++k){
      double sumjk = 0.0;
      // sum over phases j
      for(int j=0;j<F;++j){
        const auto& Mj = (j==0 ? MV : ML);
        double β = (j==0 ? betaV : betaL);
        // sum over components i,p
        for(int i=0;i<C;++i){
          if (elementMatrix[ell][i] == 0.0) continue;
          for(int p=0;p<C;++p){
            if (elementMatrix[k][p] == 0.0) continue;
            sumjk += elementMatrix[ell][i]
                  * (β * Mj[i][p])
                  * elementMatrix[k][p];
          }
        }
      }
      Acoef[ell*N + k] = sumjk;//全局矩阵中第 ell行k列的元素
    }
  } 

  // 2) 左下/右上块 A X
  // X is diag of ones per phase: A * (vector of ones) = sum over columns of A
  for(int ell = 0; ell < E; ++ell){
    // vapor 相 (第 E+0 列)
    double AX_v = 0.0;
    for(int i=0;i<C;++i){
      // A_{ell,i} * xV_i
      AX_v += elementMatrix[ell][i] * xV[i];
    }
    Acoef[ell*N + (E+0)] = AX_v;
    Acoef[(E+0)*N + ell] = AX_v;  // 对称

    // liquid 相 (第 E+1 列)
    double AX_l = 0.0;
    for(int i=0;i<C;++i){
      AX_l += elementMatrix[ell][i] * xL[i];
    }
    Acoef[ell*N + (E+1)] = AX_l;
    Acoef[(E+1)*N + ell] = AX_l;
  }
  std::cout << "betaV=" << betaV << ", betaL=" << betaL << std::endl;
  std::cout << "xV: "; for (auto x : xV) std::cout << x << " "; std::cout << std::endl;
  // 3) 下右 F×F 块全 0（按文献式 (4.25)）

  // 4) 构造 RHS u1 (元素守恒常数项)
  for (int ell = 0; ell < E; ++ell) {
    double val = 0.0;
    // ——— 蒸气相贡献 ———
    for (int i = 0; i < C; ++i) {
      double A_li = elementMatrix[ell][i];
      if (A_li == 0.0) continue;
      for (int p = 0; p < C; ++p) {
        val += A_li
             * (betaV * MV[i][p] * (muV[p] / RT));
      }
    }
    // ——— 液相贡献 ———
    for (int i = 0; i < C; ++i) {
      double A_li = elementMatrix[ell][i];
      if (A_li == 0.0) continue;
      for (int p = 0; p < C; ++p) {
        val += A_li
             * (betaL * ML[i][p] * (muL[p] / RT));
      }
    }
    rhs[ell] = val;
  }


  // 5) 构造 RHS u2 (还原自由能)
  // red_j = sum_ell ( A * x_j )_ell * (μ_j/RT)
  // but simpler: x_i = n_i/β, ∑_i x_i μ_i/RT, same as before
  double redV=0.0, redL=0.0;
  for(int i=0;i<C;++i){
    redV += (moleNumbersV[i]/betaV) * (muV[i]/RT);
    redL += (moleNumbersL[i]/betaL) * (muL[i]/RT);
  }
  rhs[E+0] = redV;
  rhs[E+1] = redL;


  // // 预计算 n_tot 和元素残差 r_elem = A*(F - (nV+nL))
  // std::vector<double> n_tot(C);
  // for (size_t i = 0; i < C; ++i) n_tot[i] = moleNumbersV[i] + moleNumbersL[i];

  // // u1: for each element ell
  // for (size_t ell = 0; ell < E; ++ell) {
  //     double r_elem = 0.0;    // A*(F - n)_ell
  //     double mu_term = 0.0;   // sum_i A_{ell,i} * [ betaV*(muV_i/RT) + betaL*(muL_i/RT) ]
  //     for (size_t i = 0; i < C; ++i) {
  //         const double A_li = elementMatrix[ell][i];
  //         r_elem += A_li * (feedComposition[i] - n_tot[i]);
  //         mu_term += A_li * (betaV * (muV[i]/RT) + betaL * (muL[i]/RT));
  //     }
  //     rhs[ell] = r_elem + mu_term;
  // }

  // // u2: each phase j uses  beta_j * (x_j · mu_j/RT)
  // double redV = 0.0, redL = 0.0;
  // for (size_t i = 0; i < C; ++i) {
  //     redV += xV[i] * (muV[i]/RT);
  //     redL += xL[i] * (muL[i]/RT);
  // }
  // rhs[E + 0] = betaV * redV;
  // rhs[E + 1] = betaL * redL;
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
  const std::vector<std::vector<double>>& MV,
  const std::vector<std::vector<double>>& ML,
  const std::vector<double>& muV,
  const std::vector<double>& muL,
  const std::vector<double>& nV,
  const std::vector<double>& nL,
  const std::vector<double>& Lambda,
  const std::vector<double>& deltaBeta,
  std::vector<double>& dnV,
  std::vector<double>& dnL) const
{
  const size_t C = nV.size();
  const size_t E = elementMatrix.size();
  (void)E; // 仅用于一致性

  dnV.assign(C, 0.0);
  dnL.assign(C, 0.0);

  const double RT = R_CONST * temperature;
  const double betaV = std::accumulate(nV.begin(), nV.end(), 0.0);
  const double betaL = std::accumulate(nL.begin(), nL.end(), 0.0);

  // A^T * Lambda
  std::vector<double> lambdaComponent(C, 0.0);
  for (size_t p = 0; p < C; ++p) {
    for (size_t ell = 0; ell < E; ++ell) {
      lambdaComponent[p] += elementMatrix[ell][p] * Lambda[ell];
    }
  }
  // diff_j = A^T λ − μ_j/RT
  std::vector<double> diffV(C), diffL(C);
  for (size_t p = 0; p < C; ++p) {
    diffV[p] = lambdaComponent[p] - muV[p] / RT;
    diffL[p] = lambdaComponent[p] - muL[p] / RT;
  }

  // 回代：Δn_i,j = x_i,j Δβ_j + β_j * M_j * diff_j
  for (size_t i = 0; i < C; ++i) {
    const double xVi = nV[i] / betaV;
    const double xLi = nL[i] / betaL;
    double combV = 0.0, combL = 0.0;
    for (size_t p = 0; p < C; ++p) {
      combV += MV[i][p] * diffV[p];
      combL += ML[i][p] * diffL[p];
    }
    dnV[i] = xVi * deltaBeta[0] + betaV * combV;
    dnL[i] = xLi * deltaBeta[1] + betaL * combL;
  }

  // 与原实现保持一致的打印
  std::cout << "\n  dnV: ";
  for (double dv : dnV) std::cout << dv << " ";
  std::cout << "\n  dnL: ";
  for (double dv : dnL) std::cout << dv << " ";
  std::cout << "\n" << std::endl;
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
  const std::vector<double>& muV,
  const std::vector<double>& muL,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& nV,
  const std::vector<double>& nL,
  const std::vector<double>& feedComposition) const
{
  RandFlash::ConvergenceInfo info{};
  int C = nV.size();
  int E = elementMatrix.size();
  // max |muV - muL|
  double max_mu_diff = 0.0;
  for (size_t i = 0; i < muV.size(); ++i) {
    max_mu_diff = std::max(max_mu_diff, std::fabs(muV[i] - muL[i]));
  }
  // 元素守恒 L∞ 范数
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

  // 与原日志一致的输出
  std::cout << "Iter "<< iternumber <<" | max_mu_diff: " << max_mu_diff << std::endl;
  std::cout << "Iter "<< iternumber <<" | element error: " << elem_err << "\n" << std::endl;
  
  // // 假设 elementMatrix 是 E×C，feed 是原始进料摩尔数
  // std::vector<double> nTot(C);
  // for (size_t i = 0; i < C; ++i) nTot[i] = nV[i] + nL[i];

  // for (size_t e = 0; e < E; ++e) {
  //     double lhs = 0.0;
  //     double rhs = 0.0;
  //     for (size_t i = 0; i < C; ++i) {
  //         lhs += elementMatrix[e][i] * nTot[i];
  //         rhs += elementMatrix[e][i] * feedComposition[i];
  //     }
  //     double err = lhs - rhs;
  //     std::cout << "element " << e << " residual = " << err << "\n";
  // }

  info.max_mu_diff = max_mu_diff;
  info.elem_error  = elem_err;
  info.converged   = (max_mu_diff < tol && elem_err < 1e-6);
  return info;
}

randflash::FlashResult RandFlash::solveTwoPhase(
  const thermo::PhaseState& state,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& initialVaporComposition,   // 可空
  const std::vector<double>& initialLiquidComposition,  // 可空
  int maxIter,
  double tol)
{
  const size_t C = state.moleNumbers.size();
  const size_t E = elementMatrix.size();

  // 计算总进料摩尔数（目前只用于 debug）
  double totalFeed = std::accumulate(state.moleNumbers.begin(),
                                     state.moleNumbers.end(), 0.0);

  // 1) 初始猜 nV, nL（优先使用提供的初始组成，否则使用平分策略）
  std::vector<double> nV(C), nL(C);
  auto init = initializeTwoPhase(
      state.moleNumbers,
      state.Temperature,
      state.Pressure,
      initialVaporComposition,
      initialLiquidComposition,
      /*margin=*/1e-3);
  nV = std::move(init.nV);
  nL = std::move(init.nL);

  std::cout << "Calculated betaV: " << init.betaV
            << "\n初始气相摩尔分率 betaV/(betaV+betaL) = "
            << (init.betaV / (init.betaV + init.betaL))
            << std::endl;

  // 2) 构造两相的 PhaseState（phaseFlag 用 backend 提供的 vapor/liquid 标志）
  thermo::PhaseState vap{state.Temperature, state.Pressure, nV, thermo_.vaporPhaseFlag()};
  thermo::PhaseState liq{state.Temperature, state.Pressure, nL, thermo_.liquidPhaseFlag()};
  // std::cout<<"构造两相PhaseState的 phaseFlag  "<<thermo_.vaporPhaseFlag()<<"  "<<thermo_.liquidPhaseFlag()<<std::endl;
  
  FlashResult result;
  result.pressure        = state.Pressure;
  result.temperature     = state.Temperature;
  result.feedComposition = state.moleNumbers;
  result.success         = false;  // 默认失败，收敛时再置 true

  // 分配临时存储（在循环外分配，循环内重复复用）
  std::vector<std::vector<double>> mV, mL;     // 局部 m_j
  std::vector<double>              muV, muL;   // 局部 μ
  std::vector<std::vector<double>> MV, ML;     // m_j 反演后的 M_j
  std::vector<double> Acoef, rhs, sol;         // 全局系统 A * x = rhs
  std::vector<double> Lambda;                  // 元素平衡拉格朗日乘子
  std::vector<double> deltaBeta(2);            // Δβ_v, Δβ_l
  std::vector<double> dnV(C), dnL(C);          // Δn^V, Δn^L

  std::cout << "进入迭代循环, maxIter = " << maxIter << std::endl;

  try {
    for (int iter = 0; iter < maxIter; ++iter) {
      // --- 2.1 更新 PhaseState 中的 n ---
      vap.moleNumbers = nV;
      liq.moleNumbers = nL;

      std::cout << "开始组装局部矩阵, C = "
                << state.moleNumbers.size()
                << " (迭代 " << (iter+1) << ")\n";

      // --- 2.2 组装两相的局部 Hessian m_j 和 μ ---
      assembleLocalJacobian(vap, mV, muV);
      assembleLocalJacobian(liq, mL, muL);
      // std::cout<<"气相组成化学势为："<<std::endl;
      // for (size_t i = 0; i < C; ++i) {
      //       std::cout << muV[i] <<"  " ;
      // }
      // std::cout<<"\n"<<"液相组成化学势为："<<std::endl;;
      // for (size_t i = 0; i < C; ++i) {
      //   std::cout << muL[i] << "  ";
      // }
      std::cout << "局部矩阵组装完成, 迭代 " << (iter+1) << std::endl;

      // === 2.3 相内 Hessian 修正 ===
      double betaV = std::accumulate(nV.begin(), nV.end(), 0.0);
      double betaL = std::accumulate(nL.begin(), nL.end(), 0.0);
      std::vector<double> xV(C), xL(C);
      for (size_t i = 0; i < C; ++i) {
        xV[i] = (betaV > 0.0) ? nV[i] / betaV : 1.0 / double(C);
        xL[i] = (betaL > 0.0) ? nL[i] / betaL : 1.0 / double(C);
      }

      PhaseFixOptions pfx;
      auto fixV = fix_phase_hessian_one_phase(linearSolver_,xV, mV, pfx);
      if (fixV.applied) {
        mV = fixV.m_fixed;
        std::cout << "[PhaseFix] Vapor: lam_min "
                  << fixV.lam_min_before << " -> "
                  << fixV.lam_min_after << "\n";
      }
      auto fixL = fix_phase_hessian_one_phase(linearSolver_,xL, mL, pfx);
      if (fixL.applied) {
        mL = fixL.m_fixed;
        std::cout << "[PhaseFix] Liquid: lam_min "
                  << fixL.lam_min_before << " -> "
                  << fixL.lam_min_after << "\n";
      }

      // 2.4 反演 m_j -> M_j
      MV = invert(mV,linearSolver_);
      ML = invert(mL,linearSolver_);
      std::cout << "局部矩阵反演完成, 迭代 " << (iter+1) << std::endl;
      // for (size_t i = 0; i < C; ++i) {
      //   for (size_t k = 0; k < C; ++k) {
      //       std::cout << MV[i][k] << "  ";  // 同一行打印
      //   }
      //   std::cout << std::endl;  // 换行
      // } 
      
      // 3) 全局系统装配
      assembleGlobalSystem(
        state.Temperature,
        MV, ML, muV, muL,
        nV, nL,
        elementMatrix,
        state.moleNumbers,   // 进料 n^F
        Acoef,
        rhs
      );
      std::cout<<"rhs: ";
      for(size_t i =0 ; i < rhs.size() ; ++i){
        std::cout<<rhs[i]<<"  ";
      }
      std::cout<<std::endl;
      std::cout << "全局系统装配完成, 迭代 " << (iter+1) << std::endl;

      // 4) 解线性系统
      double lin_resid = 0.0;
      sol = solveGlobalLinearSystem(Acoef, rhs, &lin_resid);
      std::cout << "线性系统求解完成, 迭代 " << (iter+1) << std::endl;

      // 拆解解向量：[Δλ_1…Δλ_E, Δβ_v, Δβ_l]
      Lambda.assign(sol.begin(), sol.begin() + E);
      deltaBeta[0] = sol[E];
      deltaBeta[1] = sol[E+1];

      std::cout << "  Lambda: ";
      for (double dl : Lambda) std::cout << dl << " ";
      std::cout << "\n  deltaBeta_v=" << deltaBeta[0]
                << "  deltaBeta_l=" << deltaBeta[1] << std::endl;

      // 5) 回代恢复 Δn
      {
        double betaV_dbg = std::accumulate(nV.begin(), nV.end(), 0.0);
        double betaL_dbg = std::accumulate(nL.begin(), nL.end(), 0.0);
        std::cout << "  更新前 : ";
        std::cout << "\n  nV: ";
        for (double v : nV) std::cout << v << " ";
        std::cout << "\n  nL: ";
        for (double v : nL) std::cout << v << " ";
        std::cout << "\n  betaV: " << betaV_dbg
                  << ", betaL: " << betaL_dbg;
        std::cout << "\n  xV: ";
        for (double v : nV) std::cout << v / betaV_dbg << " ";
        std::cout << "\n  xL: ";
        for (double v : nL) std::cout << v / betaL_dbg << " ";
        std::cout << std::endl;
      }

      backSubstituteDeltas(
        state.Temperature,
        elementMatrix,
        MV, ML,
        muV, muL,
        nV, nL,
        Lambda,
        deltaBeta,
        dnV, dnL
      );

      // 6) 线搜索 + 更新
      double alpha = applyUpdate(
        state.Temperature,
        muV, muL,
        dnV, dnL,
        nV, nL
      );

      // 7) 收敛判断
      auto conv = checkConvergence(
        iter + 1,
        tol,
        muV, muL,
        elementMatrix,
        nV, nL,
        state.moleNumbers
      );

      if (conv.converged) {
        double betaV_fin = std::accumulate(nV.begin(), nV.end(), 0.0);
        double betaL_fin = std::accumulate(nL.begin(), nL.end(), 0.0);
        result.success           = true;
        result.vaporFraction     = betaV_fin / (betaV_fin + betaL_fin);
        result.vaporComposition  = nV;
        result.liquidComposition = nL;
        result.iterations        = iter + 1;
        result.convergenceError  = conv.max_mu_diff;
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
  result.success = false;
  return result;
}


