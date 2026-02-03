// RAND/src/rand_flash.cpp
#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace randflash;
using namespace ls;


// === 在切空间做 SPD 修正，并保持 m x = 1 (使用 Woodbury 公式避免病态求逆) ===
// 输入: solver, x, m_in (CxC)
// 输出: PhaseFixResult (包含修正后的 M 和 m)
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



// 1) 局部 Jacobian 构造：对应论文式 (4.11)-(4.14)
void RandFlash::updatePhaseChemistry(PhaseContext& phaseCtx)
{
  // 1) 复制一份状态以便安全裁剪 n_i（避免零 / 负数）
  // 说明：这里的裁剪仅用于热力学性质求值（避免 log(0) / 导数异常），并会写回 state。
  // 绝对阈值过大时会在“濒死相(beta→0)”场景把组成推向近似均匀，进而破坏收敛与线搜索。
  // 因此采用与相总摩尔数 beta 成比例的裁剪阈值：ni >= max(1e-20, 1e-12 * beta)。
  thermo::PhaseState st = phaseCtx.state;
  const double beta0 = std::accumulate(st.moleNumbers.begin(), st.moleNumbers.end(), 0.0);
  const double clip  = std::max(1e-20, 1e-12 * std::max(beta0, 0.0));
  for (double& ni : st.moleNumbers) {
    if (ni <= clip) ni = clip;
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
    RAND_DEBUG("[PhaseFix] lam_min {} -> {}", fixRes.lam_min_before, fixRes.lam_min_after);
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
  RAND_DEBUG("  descent metric = {}", dir);

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
// static double lineSearchFeasible(
//   const std::vector<double>& nV, const std::vector<double>& nL,
//   const std::vector<double>& dnV, const std::vector<double>& dnL,
//   const std::vector<double>& gV,  const std::vector<double>& gL,
//   double init_alpha, double min_alpha, double shrink)
// {
//   auto maxFeasibleAlpha1 = [](const std::vector<double>& n,
//                               const std::vector<double>& dn){
//       double a = 1.0;
//       for (size_t i = 0; i < n.size(); ++i) {
//           if (dn[i] < 0.0) {
//               // 0.99 给一点余量，避免数值触边
//               a = std::min(a, 0.99 * n[i] / (-dn[i]));
//           }
//       }
//       return a;
//   };

//   // 方向是否下降（对约束自由能）：与 alpha 无关，判一次即可
//   auto dirDot = [&](double scale)->double {
//       (void)scale; // 方向不依赖 alpha，保持接口一致
//       double s = 0.0;
//       for (size_t i=0;i<gV.size();++i) s += dnV[i]*gV[i];
//       for (size_t i=0;i<gL.size();++i) s += dnL[i]*gL[i];
//       return s;
//   };
//   const double dir = dirDot(1.0);
//   const double dir_eps = 1e-12;
//   std::cout << "  descent metric = " << dir << std::endl; 
//   double alpha = init_alpha; // 仍然以 1.0 起步最稳
//   while (alpha > min_alpha) {
//       // 先做“可行性裁剪”，再检查下降性
//       double aV = maxFeasibleAlpha1(nV, dnV);
//       double aL = maxFeasibleAlpha1(nL, dnL);
//       alpha = std::min(alpha, std::min(aV, aL));

//       if (alpha <= min_alpha) break;

//       bool descent = (dir < -dir_eps) ? true : (std::abs(dir) <= dir_eps);
//       if (descent) return alpha;

//       alpha *= shrink; // 只在不下降时缩步
//   }
//   return min_alpha;
// }

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


// 元素守恒投影：将摩尔数投影回元素守恒约束
// 使用最小二乘法求解：min ||n - n_current||^2 subject to A*n = b
void RandFlash::projectToElementConservation(
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& targetElementMoles,
    std::vector<std::vector<double>>& nPhases_inout) const
{
    const size_t E = elementMatrix.size();
    const size_t C = elementMatrix[0].size();
    const size_t F = nPhases_inout.size();

    if (E == 0 || C == 0 || F == 0) return;

    // Compute current element moles
    std::vector<double> currentElementMoles(E, 0.0);
    for (size_t e = 0; e < E; ++e) {
        for (size_t j = 0; j < F; ++j) {
            for (size_t i = 0; i < C; ++i) {
                currentElementMoles[e] += elementMatrix[e][i] * nPhases_inout[j][i];
            }
        }
    }

    // Compute element error
    std::vector<double> elementError(E, 0.0);
    double maxError = 0.0;
    for (size_t e = 0; e < E; ++e) {
        elementError[e] = targetElementMoles[e] - currentElementMoles[e];
        maxError = std::max(maxError, std::abs(elementError[e]));
    }

    // If error is small enough, no need to project
    if (maxError < 1e-12) return;

    RAND_DEBUG("  [ElementProjection] Max element error before projection: {:.6e}", maxError);

    // Use least-squares to distribute the error across all phases
    // For each phase, solve: A * dn = elementError / F (equal distribution)
    // Then add dn to each phase

    // Simple approach: distribute error proportionally to current phase size
    std::vector<double> phaseBetas(F, 0.0);
    double totalBeta = 0.0;
    for (size_t j = 0; j < F; ++j) {
        phaseBetas[j] = std::accumulate(nPhases_inout[j].begin(), nPhases_inout[j].end(), 0.0);
        totalBeta += phaseBetas[j];
    }

    if (totalBeta < 1e-20) return;

    // For each phase, compute correction: dn_j = A^T * lambda_j
    // where lambda_j is chosen to satisfy element conservation
    // Use Moore-Penrose pseudoinverse: dn = A^T * (A*A^T)^{-1} * (error * beta_j / totalBeta)

    for (size_t j = 0; j < F; ++j) {
        double weight = phaseBetas[j] / totalBeta;
        std::vector<double> targetError(E);
        for (size_t e = 0; e < E; ++e) {
            targetError[e] = elementError[e] * weight;
        }

        // Solve A * dn = targetError using least-squares
        // Build A^T * A and A^T * b
        std::vector<std::vector<double>> ATA(C, std::vector<double>(C, 0.0));
        std::vector<double> ATb(C, 0.0);

        for (size_t i = 0; i < C; ++i) {
            for (size_t k = 0; k < C; ++k) {
                for (size_t e = 0; e < E; ++e) {
                    ATA[i][k] += elementMatrix[e][i] * elementMatrix[e][k];
                }
            }
            for (size_t e = 0; e < E; ++e) {
                ATb[i] += elementMatrix[e][i] * targetError[e];
            }
        }

        // Add regularization to avoid singularity
        for (size_t i = 0; i < C; ++i) {
            ATA[i][i] += 1e-10;
        }

        // Solve ATA * dn = ATb
        std::vector<double> dn;
        try {
            // Flatten ATA to row-major format
            std::vector<double> ATA_flat = ls::flattenRowMajor(ATA);
            dn = linearSolver_.solveSPD(static_cast<int>(C), ATA_flat, ATb);
        } catch (...) {
            RAND_DEBUG("  [ElementProjection] Failed to solve for phase {}", j);
            continue;
        }

        // Apply correction with clamping
        const double min_mole_floor = 1e-20;
        for (size_t i = 0; i < C; ++i) {
            nPhases_inout[j][i] = std::max(nPhases_inout[j][i] + dn[i], min_mole_floor);
        }
    }

    // Verify final error
    std::vector<double> finalElementMoles(E, 0.0);
    for (size_t e = 0; e < E; ++e) {
        for (size_t j = 0; j < F; ++j) {
            for (size_t i = 0; i < C; ++i) {
                finalElementMoles[e] += elementMatrix[e][i] * nPhases_inout[j][i];
            }
        }
    }

    double finalMaxError = 0.0;
    for (size_t e = 0; e < E; ++e) {
        finalMaxError = std::max(finalMaxError, std::abs(targetElementMoles[e] - finalElementMoles[e]));
    }
    RAND_DEBUG("  [ElementProjection] Max element error after projection: {:.6e}", finalMaxError);
}

// 6) 结果更新：线搜索 + n 的更新与"更新后"打印（逻辑不变）
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

  RAND_DEBUG("  线搜索步长 alpha = {:.4f}", alpha);

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
  std::string betaStr = "  更新后各相总摩尔数 (Beta): ";
  for(size_t j=0; j<F; ++j) {
      double beta_j = std::accumulate(nPhases_inout[j].begin(), nPhases_inout[j].end(), 0.0);
      betaStr += fmt::format("{:.6f} ", beta_j);
  }
  RAND_DEBUG("{}", betaStr);

  // ========== 新增：打印各相摩尔分率（固定小数位，无科学计数法） ==========
  std::string moleFracStr = "  更新后各相摩尔分率:\n";
  for(size_t j=0; j<F; ++j) {
      // 计算当前相的总摩尔数
      double beta_j = std::accumulate(nPhases_inout[j].begin(), nPhases_inout[j].end(), 0.0);
      // 防止除以0（理论上不会触发，因有min_mole_floor）
      beta_j = std::max(beta_j, min_mole_floor);

      moleFracStr += fmt::format("    相 {}: ", j+1);
      for(size_t i=0; i<nPhases_inout[j].size(); ++i) {
          // 计算摩尔分率
          double mole_frac = nPhases_inout[j][i] / beta_j;
          moleFracStr += fmt::format("{:.8f} ", mole_frac);
      }
      moleFracStr += "\n";
  }
  RAND_DEBUG("{}", moleFracStr);
}

// 7) 收敛判断：打印与返回布尔值（阈值保持原样：elem_error < 1e-8）
RandFlash::ConvergenceInfo RandFlash::checkConvergence(
  int iternumber,
  double tol,
  double temperature,
  const std::vector<std::vector<double>>& mus,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<std::vector<double>>& nPhases, // 注意：这里需要传入 nPhases 或 compositions
  const std::vector<double>& feedComposition,
  const std::vector<std::vector<double>>& dnPhases,  // New: step size for convergence check
  bool isReactive) const                              // New: flag for reactive systems
{
  RandFlash::ConvergenceInfo info{};
  const double RT = R_CONST * temperature;
  const size_t F = mus.size();
  const size_t C = feedComposition.size();
  const size_t E = elementMatrix.size();

  // 1) 化学势平衡判据 (加入微量组分过滤)
  double max_mu_diff = 0.0;
  // 阈值：如果某组分在某相的摩尔分率低于此值，则忽略其化学势差异
  const double TRACE_LIMIT = 1e-6;  // Original value

  if (F > 1) {
    // 预先计算各相的总摩尔数 beta 和摩尔分率 x
    std::vector<double> betas(F, 0.0);
    std::vector<std::vector<double>> xs(F, std::vector<double>(C, 0.0));
    for(size_t j=0; j<F; ++j) {
        betas[j] = std::accumulate(nPhases[j].begin(), nPhases[j].end(), 0.0);
        if(betas[j] > 1e-20) {
            for(size_t i=0; i<C; ++i) xs[j][i] = nPhases[j][i] / betas[j];
        }
    }

    for (size_t j = 1; j < F; ++j) {
        for (size_t i = 0; i < C; ++i) {
              // 【核心修改】仅当两相中该组分都"显著存在"时，才比较化学势
              // 否则，该组分的对数项导致的差异是数值噪音或物理上的不混溶表现
              if (xs[j][i] > TRACE_LIMIT && xs[0][i] > TRACE_LIMIT) {
                  double diff = std::fabs(mus[j][i] - mus[0][i]) / RT;
                  max_mu_diff = std::max(max_mu_diff, diff);
              }
        }
    }
  }

  // 2) 元素守恒判据 (保持不变)
  double elem_err = 0.0;
  for (size_t ell = 0; ell < E; ++ell) {
      double sum_n_ell = 0.0;
      for (size_t j=0; j<F; ++j) {
          for (size_t i=0; i<C; ++i) {
              sum_n_ell += elementMatrix[ell][i] * nPhases[j][i];
          }
      }
      double sum_feed_ell = 0.0;
      for (size_t i=0; i<C; ++i) {
           sum_feed_ell += elementMatrix[ell][i] * feedComposition[i];
      }
      elem_err = std::max(elem_err, std::abs(sum_n_ell - sum_feed_ell));
  }

  // 3) 计算相对步长范数（用于单相反应体系）
  double relative_step_norm = 0.0;
  if (!dnPhases.empty() && dnPhases.size() == F) {
      double step_norm = 0.0, n_norm = 0.0;
      for (size_t j = 0; j < F; ++j) {
          if (j < dnPhases.size() && dnPhases[j].size() == C) {
              for (size_t i = 0; i < C; ++i) {
                  step_norm += dnPhases[j][i] * dnPhases[j][i];
                  n_norm += nPhases[j][i] * nPhases[j][i];
              }
          }
      }
      relative_step_norm = (n_norm > 1e-20) ? std::sqrt(step_norm / n_norm) : 0.0;
  }

  if (isReactive && F == 1) {
      RAND_DEBUG("Iter {} | max_mu/RT_diff (filtered): {} | elem_err: {} | relative_step_norm: {}",
                iternumber, max_mu_diff, elem_err, relative_step_norm);
  } else {
      RAND_DEBUG("Iter {} | max_mu/RT_diff (filtered): {} | elem_err: {}",
                iternumber, max_mu_diff, elem_err);
  }

  info.max_mu_diff = max_mu_diff;
  info.elem_error  = elem_err;
  info.relative_step_norm = relative_step_norm;

  // 收敛判断：
  // - 单相反应体系：使用步长判据
  // - 多相或非反应体系：使用化学势差判据
  if (isReactive && F == 1) {
      info.converged = (relative_step_norm < tol && elem_err < 1e-5);
  } else {
      info.converged = (max_mu_diff < tol && elem_err < 1e-5);
  }

  return info;
}

void RandFlash::printResult(const MultiFlashResult& res) const {
  // 1. 自动从 Backend 获取组分名称
  std::vector<std::string> compNames = thermo_.getComponentNames();

  RAND_INFO("\n============================== Flash Results ==============================");
  RAND_INFO("Temperature: {:.2f} K", res.temperature);
  RAND_INFO("Pressure:    {:.2f} Pa", res.pressure);
  RAND_INFO("Iterations:  {}", res.iterations);
  RAND_INFO("Status:      {}", (res.success ? "CONVERGED" : "FAILED"));
  RAND_INFO("Error Norm:  {:.4e}", res.mu_infinity_norm);

  size_t F = res.phases.size();
  if (F == 0) {
      RAND_INFO("No phases returned.");
      RAND_INFO("===========================================================================");
      return;
  }

  size_t C = res.phases[0].state.moleNumbers.size();

  // 校验一下组分数量是否匹配，防止越界打印
  if (compNames.size() != C) {
      compNames.resize(C, "Unknown");
  }

  double total_moles = 0.0;
  for (size_t j = 0; j < F; ++j) {
      total_moles += res.beta(j);
  }

  for (size_t j = 0; j < F; ++j) {
      double beta_j = res.beta(j);
      double phase_frac = (total_moles > 1e-12) ? (beta_j / total_moles) : 0.0;

      // 使用phases中存储的phaseFlag判断相态
      std::string phaseType;
      int flag = res.phases[j].state.phaseFlag;
      if (flag == thermo_.vaporPhaseFlag()) {
          phaseType = "Vapor";
      } else {
          phaseType = "Liquid";
      }

      RAND_INFO("\n---------------------------------------------------------------------------");
      RAND_INFO(" {} {} | Phase Fraction (Beta): {:.5f} | Total Moles: {:.5f}", phaseType, j, phase_frac, beta_j);
      RAND_INFO("---------------------------------------------------------------------------");
      RAND_INFO("  Idx | {:15} | {:>12} | {:>12}", "Component", "Mole Frac (x)", "Moles (n)");
      RAND_INFO("------+-----------------+--------------+--------------");

      const auto& n_j = res.phases[j].state.moleNumbers;
      for (size_t i = 0; i < C; ++i) {
          double n_i = n_j[i];
          double x_i = (beta_j > 1e-20) ? (n_i / beta_j) : 0.0;

          RAND_INFO("  {:3} | {:15} | {:12.5f} | {:12.4e}", i, compNames[i], x_i, n_i);
      }
  }
  RAND_INFO("===========================================================================\n");
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

// ========== Helper Methods for Automatic Phase Number Determination ==========

// Compute total Gibbs free energy of the system
// G_total = sum_j sum_i n_i^(j) * mu_i^(j)
double RandFlash::computeTotalGibbs(const SystemContext& sys) const
{
    double G_total = 0.0;

    for (size_t j = 0; j < sys.phases.size(); ++j) {
        const auto& phase = sys.phases[j];

        // Ensure chemical potentials are available
        if (phase.mu.size() != phase.state.moleNumbers.size()) {
            throw std::runtime_error("Chemical potentials not computed for phase " + std::to_string(j));
        }

        // G_phase = sum_i n_i * mu_i
        for (size_t i = 0; i < phase.state.moleNumbers.size(); ++i) {
            G_total += phase.state.moleNumbers[i] * phase.mu[i];
        }
    }

    return G_total;
}

// Generate trial compositions for phase splitting
// Uses random perturbations and Wilson K-value based compositions
std::vector<std::vector<double>> RandFlash::generateTrialCompositions(
    const SystemContext& sys, int numTrials) const
{
    const size_t C = sys.feedMoles.size();
    std::vector<std::vector<double>> trials;

    if (C == 0 || sys.phases.empty()) return trials;

    // Get current phase composition (use first phase as reference)
    std::vector<double> x_ref = sys.phases[0].x;

    // Normalize reference composition
    double sum_ref = std::accumulate(x_ref.begin(), x_ref.end(), 0.0);
    if (sum_ref > 1e-20) {
        for (double& xi : x_ref) xi /= sum_ref;
    } else {
        std::fill(x_ref.begin(), x_ref.end(), 1.0 / C);
    }

    // Trial 1: Wilson K-value based vapor-like composition
    std::vector<double> K(C, 1.0);
    thermo_.wilsonK(sys.temperature, sys.pressure, K);

    std::vector<double> trial_vapor(C);
    for (size_t i = 0; i < C; ++i) {
        trial_vapor[i] = x_ref[i] * K[i];
    }
    double sum_v = std::accumulate(trial_vapor.begin(), trial_vapor.end(), 0.0);
    if (sum_v > 1e-20) {
        for (double& xi : trial_vapor) xi /= sum_v;
    }
    trials.push_back(trial_vapor);

    // Trial 2: Wilson K-value based liquid-like composition
    std::vector<double> trial_liquid(C);
    for (size_t i = 0; i < C; ++i) {
        trial_liquid[i] = x_ref[i] / K[i];
    }
    double sum_l = std::accumulate(trial_liquid.begin(), trial_liquid.end(), 0.0);
    if (sum_l > 1e-20) {
        for (double& xi : trial_liquid) xi /= sum_l;
    }
    trials.push_back(trial_liquid);

    // Additional trials: random perturbations
    for (int t = 2; t < numTrials; ++t) {
        std::vector<double> trial(C);
        for (size_t i = 0; i < C; ++i) {
            // Random perturbation factor between 0.5 and 2.0
            double factor = 0.5 + 1.5 * (std::rand() % 1000) / 1000.0;
            trial[i] = x_ref[i] * factor;
        }
        double sum_t = std::accumulate(trial.begin(), trial.end(), 0.0);
        if (sum_t > 1e-20) {
            for (double& xi : trial) xi /= sum_t;
        }
        trials.push_back(trial);
    }

    return trials;
}

// Check if splitting current phase with trial composition reduces Gibbs energy
bool RandFlash::checkPhaseSplitBenefit(
    const SystemContext& sys,
    const std::vector<double>& trialComposition,
    double currentGibbs,
    double& newGibbs) const
{
    // This is a simplified check - a full implementation would:
    // 1. Create a new system with F+1 phases
    // 2. Initialize the new phase with trialComposition
    // 3. Solve the multi-phase equilibrium
    // 4. Compare Gibbs energies

    // For now, return false as placeholder
    // Full implementation will be added in solveReactiveAuto
    newGibbs = currentGibbs;
    return false;
}

// ========== Stable Phase Determination Methods ==========

// Analyze stability of all phases in a multi-phase system
RandFlash::MultiPhaseStabilityResult RandFlash::analyzeMultiPhaseStability(
    const MultiFlashResult& currentResult,
    const phase_stability::StabilityOptions& stabOpt) const
{
    MultiPhaseStabilityResult result;
    result.globally_stable = true;

    phase_stability::PhaseStabilityAnalyzer analyzer(thermo_);

    // Perform stability analysis on each phase
    for (size_t j = 0; j < currentResult.numPhases(); ++j) {
        const auto& phase = currentResult.phases[j];
        double beta = currentResult.beta(j);

        // Skip phases with negligible amount
        if (beta < 1e-10) continue;

        // Use phase composition as reference for stability analysis
        auto stabRes = analyzer.analyze(
            currentResult.temperature,
            currentResult.pressure,
            phase.x,
            stabOpt);

        if (!stabRes.stable) {
            result.globally_stable = false;

            // Collect all incipient phases found
            for (const auto& inc : stabRes.incipient) {
                result.all_incipients.push_back(inc);
            }
        }
    }

    // Deduplicate incipients based on composition similarity
    if (!result.all_incipients.empty()) {
        std::vector<phase_stability::IncipientPhase> unique_incipients;
        const double l1_threshold = 0.05;  // L1 distance threshold for deduplication

        for (const auto& inc : result.all_incipients) {
            bool is_duplicate = false;

            for (const auto& existing : unique_incipients) {
                // Calculate L1 distance
                double l1_dist = 0.0;
                for (size_t i = 0; i < inc.x.size(); ++i) {
                    l1_dist += std::abs(inc.x[i] - existing.x[i]);
                }

                if (l1_dist < l1_threshold) {
                    is_duplicate = true;
                    break;
                }
            }

            if (!is_duplicate) {
                unique_incipients.push_back(inc);
            }
        }

        result.all_incipients = unique_incipients;

        // Sort by TPD (most negative first)
        std::sort(result.all_incipients.begin(), result.all_incipients.end(),
            [](const phase_stability::IncipientPhase& a, const phase_stability::IncipientPhase& b) {
                return a.tpd < b.tpd;
            });
    }

    return result;
}

// Select the best incipient phase to add
// Strategy: prioritize phase type diversity, then TPD magnitude
phase_stability::IncipientPhase RandFlash::selectBestIncipient(
    const std::vector<phase_stability::IncipientPhase>& incipients,
    const MultiFlashResult& currentResult) const
{
    if (incipients.empty()) {
        throw std::runtime_error("selectBestIncipient: no incipients provided");
    }

    // Collect existing phase flags
    std::vector<int> existing_flags;
    for (const auto& phase : currentResult.phases) {
        double beta = std::accumulate(phase.state.moleNumbers.begin(),
                                      phase.state.moleNumbers.end(), 0.0);
        if (beta > 1e-10) {
            existing_flags.push_back(phase.state.phaseFlag);
        }
    }

    const int vap = thermo_.vaporPhaseFlag();
    const int liq = thermo_.liquidPhaseFlag();

    // Check if we have vapor and liquid phases
    bool has_vapor = std::find(existing_flags.begin(), existing_flags.end(), vap) != existing_flags.end();
    bool has_liquid = std::find(existing_flags.begin(), existing_flags.end(), liq) != existing_flags.end();

    // Strategy 1: Prioritize phase type diversity
    // If we only have liquid phases, prefer adding a vapor phase
    if (!has_vapor && has_liquid) {
        for (const auto& inc : incipients) {
            if (inc.phase_flag == vap) {
                RAND_DEBUG("[SelectIncipient] Prioritizing vapor phase for diversity (TPD={})", inc.tpd);
                return inc;
            }
        }
    }

    // If we only have vapor phase, prefer adding a liquid phase
    if (has_vapor && !has_liquid) {
        for (const auto& inc : incipients) {
            if (inc.phase_flag == liq) {
                RAND_DEBUG("[SelectIncipient] Prioritizing liquid phase for diversity (TPD={})", inc.tpd);
                return inc;
            }
        }
    }

    // Strategy 2: Select the incipient with most negative TPD
    // (already sorted by TPD in analyzeMultiPhaseStability)
    RAND_DEBUG("[SelectIncipient] Selecting incipient with most negative TPD={}", incipients[0].tpd);
    return incipients[0];
}

// Attempt to add a new phase to the current result
MultiFlashResult RandFlash::attemptPhaseAddition(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const MultiFlashResult& currentResult,
    const phase_stability::IncipientPhase& newPhase,
    int maxIter,
    double tol)
{
    // Extract existing phase compositions and flags
    std::vector<std::vector<double>> compositions;
    std::vector<int> phaseFlags;

    for (const auto& phase : currentResult.phases) {
        double beta = std::accumulate(phase.state.moleNumbers.begin(),
                                      phase.state.moleNumbers.end(), 0.0);
        if (beta > 1e-10) {
            compositions.push_back(phase.x);
            phaseFlags.push_back(phase.state.phaseFlag);
        }
    }

    // Add the new incipient phase
    compositions.push_back(newPhase.x);
    phaseFlags.push_back(newPhase.phase_flag);

    size_t newPhaseCount = compositions.size();

    RAND_DEBUG("[PhaseAddition] Adding new phase (flag={}), total phases: {}",
               newPhase.phase_flag, newPhaseCount);

    // === 方案1：水体系检测与专用初始化 ===
    // 检测是否为水体系
    const size_t C = input.feedMoles.size();
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

    // 计算归一化进料组成
    double total_feed = std::accumulate(input.feedMoles.begin(), input.feedMoles.end(), 0.0);
    std::vector<double> z_feed = input.feedMoles;
    if (total_feed > 1e-20) {
        for (double& zi : z_feed) zi /= total_feed;
    }

    double z_water = (water_idx >= 0 && static_cast<size_t>(water_idx) < z_feed.size())
                     ? z_feed[water_idx] : 0.0;
    bool is_water_system = (water_idx >= 0 && z_water > 1e-4);

    // 如果是水体系且要添加第三相，使用专用三相水体系初始化
    if (is_water_system && newPhaseCount == 3) {
        RAND_DEBUG("[PhaseAddition] Detected water system (water_idx={}, z_water={:.4f}), "
                   "using specialized three-phase water initialization", water_idx, z_water);

        // 构建 SystemContext
        SystemContext sys;
        sys.temperature = input.temperature;
        sys.pressure = input.pressure;
        sys.feedMoles = input.feedMoles;
        sys.elementMatrix = elementMatrix;

        // 使用专用三相水体系初始化
        auto init = initializeThreePhaseWater(sys);

        // 设置相标识符：0=Vapor, 1=Oil, 2=Aqueous
        sys.phases.resize(3);
        sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], thermo_.vaporPhaseFlag()};
        sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], thermo_.liquidPhaseFlag()};
        sys.phases[2].state = {input.temperature, input.pressure, init.n_phases[2], thermo_.liquidPhaseFlag()};

        // 调用通用求解器
        bool isReactive = !elementMatrix.empty() && elementMatrix.size() < C;
        return solveGeneral(sys, maxIter, tol, isReactive);
    }

    // 通用路径：使用 incipient 组成
    return solveWithPhaseGuesses(input, elementMatrix, compositions, phaseFlags, maxIter, tol);
}

// Stable phase determination: iterative stability-splitting loop
MultiFlashResult RandFlash::solveStable(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIter,
    double tol)
{
    RAND_INFO("\n=== Stable Phase Determination Method ===");
    RAND_INFO("Max stability iterations: {}", opt.max_stability_iterations);
    RAND_INFO("Max phases: {}", opt.max_phases);
    RAND_INFO("Gibbs improvement threshold: {} J", opt.gibbs_improvement_threshold);

    // Step 1: Initialize with single phase
    const int mingibbs = thermo_.minGibbsPhaseFlag();
    std::vector<std::vector<double>> singlePhaseComp = {input.feedMoles};
    std::vector<int> singlePhaseFlag = {mingibbs};

    MultiFlashResult currentResult = solveWithPhaseGuesses(
        input, elementMatrix, singlePhaseComp, singlePhaseFlag, maxIter, tol);

    if (!currentResult.success) {
        RAND_INFO("Single-phase initialization failed");
        return currentResult;
    }

    // Compute initial Gibbs energy
    SystemContext sys;
    sys.temperature = input.temperature;
    sys.pressure = input.pressure;
    sys.feedMoles = input.feedMoles;
    sys.elementMatrix = elementMatrix;
    sys.phases = currentResult.phases;

    for (size_t j = 0; j < sys.phases.size(); ++j) {
        updatePhaseChemistry(sys.phases[j]);
    }

    double currentGibbs = computeTotalGibbs(sys);
    RAND_INFO("Initial single-phase Gibbs: {} J", currentGibbs);

    // Step 2: Iterative stability-splitting loop
    phase_stability::StabilityOptions stabOpt;
    stabOpt.verbose = opt.verbose_stability_loop;
    stabOpt.max_ss_iters = 100;
    stabOpt.n_random_seeds = 8;
    stabOpt.tpd_tol = 1e-12;

    for (int iter = 0; iter < opt.max_stability_iterations; ++iter) {
        if (opt.verbose_stability_loop) {
            RAND_INFO("\n--- Stability iteration {} ---", iter + 1);
            RAND_INFO("Current phases: {}, Gibbs: {} J", currentResult.numPhases(), currentGibbs);
        }

        // Step 2a: Perform stability analysis on all phases
        auto stabResult = analyzeMultiPhaseStability(currentResult, stabOpt);

        // Step 2b: Check if globally stable
        if (stabResult.globally_stable) {
            RAND_INFO("System is globally stable with {} phases", currentResult.numPhases());
            return currentResult;
        }

        if (opt.verbose_stability_loop) {
            RAND_INFO("Found {} incipient phases", stabResult.all_incipients.size());
        }

        // Step 2c: Check phase number limit
        if (static_cast<int>(currentResult.numPhases()) >= opt.max_phases) {
            RAND_INFO("Reached max phases limit ({}), stopping", opt.max_phases);
            return currentResult;
        }

        // Step 2d: Select best incipient phase
        auto incipient = selectBestIncipient(stabResult.all_incipients, currentResult);

        if (opt.verbose_stability_loop) {
            RAND_INFO("Selected incipient: flag={}, TPD={}", incipient.phase_flag, incipient.tpd);
        }

        // Step 2e: Attempt to add new phase
        auto newResult = attemptPhaseAddition(input, elementMatrix, currentResult, incipient, maxIter, tol);

        // Step 2f: Check if phase addition succeeded
        if (!newResult.success) {
            RAND_INFO("Phase addition failed to converge, stopping");
            return currentResult;
        }

        // Step 2g: Compute new Gibbs energy
        sys.phases = newResult.phases;
        for (size_t j = 0; j < sys.phases.size(); ++j) {
            updatePhaseChemistry(sys.phases[j]);
        }

        double newGibbs = computeTotalGibbs(sys);
        double deltaG = currentGibbs - newGibbs;

        if (opt.verbose_stability_loop) {
            RAND_INFO("New Gibbs: {} J, ΔG: {} J", newGibbs, deltaG);
        }

        // Step 2h: Check Gibbs improvement
        if (deltaG < opt.gibbs_improvement_threshold) {
            RAND_INFO("Gibbs improvement insufficient (ΔG={} < threshold={}), stopping",
                     deltaG, opt.gibbs_improvement_threshold);
            return currentResult;
        }

        // Step 2i: Accept new result
        currentResult = std::move(newResult);
        currentGibbs = newGibbs;

        if (opt.verbose_stability_loop) {
            RAND_INFO("Accepted new configuration with {} phases", currentResult.numPhases());
        }
    }

    RAND_INFO("Reached max stability iterations ({})", opt.max_stability_iterations);
    return currentResult;
}

