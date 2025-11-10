#include "rand_flash.hpp"
#include <numeric>
#include <cmath>
using namespace randflash;
using namespace ls;

// 保证逐分量守恒 + 组分单纯形精确成立 (暂无使用)
inline void enforce_component_balance_and_simplex(
  Eigen::VectorXd& nV, Eigen::VectorXd& nL,
  const Eigen::VectorXd& nFeed, // 全流程不变的进料逐分量摩尔数
  double* betaV_out = nullptr, double* betaL_out = nullptr,
  Eigen::VectorXd* xV_out = nullptr, Eigen::VectorXd* xL_out = nullptr)
{
  const int nc = static_cast<int>(nV.size());
  assert(nc == nL.size() && nc == nFeed.size());

  // 1) 逐分量裁剪，避免数值微负/微超
  for (int i = 0; i < nc; ++i) {
      // 把 nV 限制在 [0, nF]
      if (nV[i] < 0.0) nV[i] = 0.0;
      if (nV[i] > nFeed[i]) nV[i] = nFeed[i];
      // 由守恒强制 nL
      nL[i] = nFeed[i] - nV[i];
      // 再保一次非负（理论上不会触发，防卫式）
      if (nL[i] < 0.0) { nV[i] += nL[i]; nL[i] = 0.0; } // 把负数挪到另一相
      if (nV[i] < 0.0) { nL[i] += nV[i]; nV[i] = 0.0; }
  }

  // 2) 精确计算 beta，并构造单纯形上的 x
  const double betaV = nV.sum();
  const double betaL = nL.sum();

  if (betaV_out) *betaV_out = betaV;
  if (betaL_out) *betaL_out = betaL;

  if (xV_out) {
    *xV_out = (betaV > 0.0 
        ? (nV / betaV).eval()  // 强制除法表达式求值为VectorXd
        : Eigen::VectorXd::Zero(nc)
    );
  }
  if (xL_out) {
      *xL_out = (betaL > 0.0 
          ? (nL / betaL).eval()  // 同理，强制求值
          : Eigen::VectorXd::Zero(nc)
      );
  }
}

// === 构造切空间正交基 B: 1^T y = 0 ===
// 用 (e_i - e_C) 做初基，再 Gram–Schmidt 正交化
static Eigen::MatrixXd tangentBasis(int C){
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(C, C-1);
  for (int k=0;k<C-1;++k){ // 列 k
    B(k,   k) =  1.0;
    B(C-1,k) = -1.0;
  }
  // Gram-Schmidt
  for (int j=0;j<C-1;++j){
    for (int i=0;i<j;++i){
      double proj = B.col(i).dot(B.col(j));
      B.col(j) -= proj * B.col(i);
    }
    double nrm = B.col(j).norm();
    if (nrm < 1e-14) { // 退化保护
      // 随机扰动再正交
      B.col(j).setRandom();
      for (int i=0;i<j;++i){
        double proj = B.col(i).dot(B.col(j));
        B.col(j) -= proj * B.col(i);
      }
      nrm = B.col(j).norm();
    }
    B.col(j) /= nrm;
  }
  // 检查 1^T B = 0
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(C);
  Eigen::RowVectorXd check = ones.transpose()*B;
  // 可加断言(略)
  return B;
}

// === 在切空间做 SPD 修正，并保持 m x = 1 ===

static PhaseFixResult fix_phase_hessian_one_phase(
  const std::vector<double>& x_std,                 // 相内 x
  const std::vector<std::vector<double>>& m_in_std, // 相内 m
  const PhaseFixOptions& opt = {})
  {
  const int C = (int)x_std.size();
  Eigen::VectorXd x = toEig(x_std);
  Eigen::MatrixXd m_in = toEig(m_in_std);
  // 对称部
  Eigen::MatrixXd ms = 0.5*(m_in + m_in.transpose());
  // 切空间特征值
  Eigen::MatrixXd B = tangentBasis(C);
  Eigen::MatrixXd mt = B.transpose() * ms * B;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es0(mt);
  double lam_min = es0.eigenvalues()(0); // 升序
  PhaseFixResult out;
  out.lam_min_before = lam_min;

  double add = 0.0;
  Eigen::MatrixXd mt_fixed = mt;
  if (lam_min < opt.eig_floor){
    add = opt.eig_floor - lam_min + opt.eps_shift;
    mt_fixed += add * Eigen::MatrixXd::Identity(C-1, C-1);
  }
  if (lam_min >= opt.eig_floor) {
    PhaseFixResult out;
    out.applied = false;
    out.lam_min_before = lam_min;
    out.lam_min_after  = lam_min;
    out.m_fixed = m_in_std;
    // 这里不必求逆，沿用你外层的 invert() 即可；如需也可在此求一下
    return out;
  }
  // 只回写切空间修正：不会动到 1 方向的块
  Eigen::MatrixXd m_tan = B * mt_fixed * B.transpose();

  // === 新：构造 K = 1 z^T + z 1^T + gamma 11^T，强制 m x = 1，且不改切空间块 ===
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(C);

  // 目标残差 r = 1 - m_tan * x
  Eigen::VectorXd r = ones - m_tan * x;

  // 分解 r = r_perp + rho * 1
  double rho = r.sum() / static_cast<double>(C);
  Eigen::VectorXd r_perp = r - rho * ones;

  // 设 z = r_perp, gamma = rho - r_perp^T x
  Eigen::VectorXd z = r_perp;
  double gamma = rho - r_perp.dot(x);

  // 组装对称修正 K
  Eigen::MatrixXd K = ones * z.transpose() + z * ones.transpose()
                    + gamma * (ones * ones.transpose());

  // 最终修正后的 m
  Eigen::MatrixXd m_fix = m_tan + K;

  // （可选）数值余量：轻微对称化，避免舍入
  m_fix = 0.5 * (m_fix + m_fix.transpose());

  // 现在检查切空间最小特征值（不会变，因为 B^T K B = 0）
  Eigen::MatrixXd mt_chk = B.transpose() * m_fix * B;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_chk(mt_chk);
  out.lam_min_after = es_chk.eigenvalues()(0);

  // 求逆
  Eigen::LDLT<Eigen::MatrixXd> ldlt(m_fix);
  if (ldlt.info()!=Eigen::Success) {
    m_fix += opt.eps_shift * Eigen::MatrixXd::Identity(C,C);
    ldlt.compute(m_fix);
  }
  Eigen::MatrixXd M_fix = ldlt.solve(Eigen::MatrixXd::Identity(C,C));

  // 不变量核验：M*1 是否等于 x（理论上应当精确成立）
  Eigen::VectorXd diff = M_fix * ones - x;
  // 允许很小数值误差
  if (diff.norm() > 1e-9) {
    // 极小对称化补偿（理论上用不到）
    M_fix -= 0.5 * (diff * ones.transpose() + ones * diff.transpose()) / static_cast<double>(C);
  }

  // 输出
  out.applied = (add > 0.0);     // 只有切空间真做了移位才算“应用”
  out.m_fixed = toStd(m_fix);
  out.M_fixed = toStd(M_fix);
  return out;
}

RandFlash::RandFlash(
    PropertyPackageType packageType,
    std::shared_ptr<material_object::Cluster> componentCluster,
    ls::LinearSolverInterface& linearSolver)
  : vaporModel_(packageType, componentCluster),
    liquidModel_(packageType, componentCluster),
    linearSolver_(linearSolver)
{
  // 如果有额外初始化,可以放在这里
}

// ---- in rand_flash.cpp ----
InitResult RandFlash::initializeTwoPhase(
  double /*pressure*/,
  double /*temperature*/,
  const std::vector<double>& feed,
  const std::vector<std::vector<double>>& /*elementMatrix*/,
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

  // 案例 A：只给了气相组成 y
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

  // 案例 B：只给了液相组成 x
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

  // 案例 C：同时给了 y 与 x —— 用最小二乘估计“共识” beta，再恢复 nV/nL
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

  // 案例 D：都没给 —— 回退到 50/50 平分
  out.betaV = 0.5*Ftot;
  out.betaL = 0.5*Ftot;
  out.nV = out.nL = feed;
  for (double &x : out.nV) x *= 0.5;
  for (double &x : out.nL) x *= 0.5;
  out.y0.resize(C); out.x0.resize(C);
  for (size_t i=0;i<C;++i) { out.y0[i] = out.nV[i]/out.betaV; out.x0[i] = out.nL[i]/out.betaL; }
  return out;
}

// 1) 局部 Jacobian 构造：对应论文式 (4.11)-(4.14)
//    m_j[i][k] = β_j * ( (1/RT) * ∂μ_i/∂n_k + 1 )
//    mu_j 已经由 chemicalPotentials 给出
void RandFlash::assembleLocalJacobian(
  double temperature,
  double pressure,
  std::vector<double>& moleNumbers,
  thermo::PropertyPackageAdapter& model,
  std::vector<std::vector<double>>& m,
  std::vector<double>& mu)
{
  for (double &ni : moleNumbers) if (ni <= 1e-10) ni = 1e-10; 

  // 1. 构造 PhaseState
  thermo::PhaseState state{ temperature, pressure, moleNumbers };

  // 2. 计算 μ 和 ∂μ/∂n
  mu = model.chemicalPotentials(state);

  auto dmun = model.dMu_dN(state);
  // for (size_t i = 0; i < moleNumbers.size(); ++i) {
  //   for (size_t k = 0; k < moleNumbers.size(); ++k) {
  //       std::cout << dmun[i][k] << "  ";  // 同一行打印
  //   }
  //   std::cout << std::endl;  // 换行
  // }    
  // 3. 构造 m_j 矩阵：m_j[i][k] = β_j*(dμ/dn/RT +1),此处 β_j = totalMoles
  double totalMoles = std::accumulate(moleNumbers.begin(), moleNumbers.end(), 0.0);
  double RT = R_CONST * temperature;
  size_t C = moleNumbers.size();
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
  double dmean = 0.0; for (size_t i=0;i<C;++i) dmean += std::abs(m[i][i]); dmean = std::max(dmean/C, 1.0);
  for (size_t i=0;i<C;++i) m[i][i] += 1e-10 * dmean;
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
  assert(static_cast<int>(Acoef.size()) == N*N);

  // Map 1D Acoef 到 Eigen 行主序矩阵
  Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
      A_mat(Acoef.data(), N, N);
  Eigen::Map<const Eigen::VectorXd> b_vec(rhs.data(), N);

  // 与原代码保持一致：使用 SVD 最小范数解
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      A_mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::VectorXd x_vec = svd.solve(b_vec);

  // 残差检测（保持原有阈值与日志）
  const double res_norm = (A_mat * x_vec - b_vec).norm();
  if (residual_out) *residual_out = res_norm;
  if (res_norm > 1e-6) {
    std::cerr << "[RandFlash] Warning: high residual = " << res_norm << "\n";
  }

  // 写回 std::vector
  std::vector<double> sol;
  sol.assign(x_vec.data(), x_vec.data() + N);
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

  info.max_mu_diff = max_mu_diff;
  info.elem_error  = elem_err;
  info.converged   = (max_mu_diff < tol && elem_err < 1e-6);
  return info;
}

FlashResult RandFlash::solveTwoPhase(
  double pressure,
  double temperature,
  const std::vector<double>& feedComposition,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<double>& initialVaporComposition,  // 新增参数
  const std::vector<double>& initialLiquidComposition, // 新增参数
  int maxIter,
  double tol)
{
  size_t C = feedComposition.size();
  size_t E = elementMatrix.size();

  // 计算总进料摩尔数
  double totalFeed = std::accumulate(feedComposition.begin(), feedComposition.end(), 0.0);
  
  // 初始猜 nV, nL：优先使用提供的初始组成,否则使用平分策略
  std::vector<double> nV(C), nL(C);
  auto init = initializeTwoPhase(
    pressure, temperature,
    feedComposition, elementMatrix,
    initialVaporComposition, initialLiquidComposition,
    /*margin=*/1e-3);
  nV = std::move(init.nV);
  nL = std::move(init.nL);

  std::cout << "Calculated beta: " << init.betaV
            << "\n初始气相摩尔分率为： " << (init.betaV / (init.betaV+init.betaL))
            << std::endl;

  FlashResult result;
  result.pressure       = pressure;
  result.temperature    = temperature;
  result.feedComposition= feedComposition;

  // 分配临时存储
  std::vector<std::vector<double>> mV, mL;     // 块 Hessian m_j
  std::vector<double> muV, muL;                // μ 向量
  std::vector<std::vector<double>> MV, ML;     // 反演后 M_j
  std::vector<double> Acoef, rhs, sol;         // 全局系统
  std::vector<double> Lambda, deltaBeta(2);
  std::vector<double> dnV(C), dnL(C);

  std::cout << "进入迭代循环,maxIter=" << maxIter << std::endl;
  try{
    for (int iter = 0; iter < maxIter; ++iter) {
      // 1) 本地 Jacobian 和 μ
      vaporModel_.setPhaseIndex(0);
      liquidModel_.setPhaseIndex(1);

      std::cout << "开始组装局部矩阵,moleNumbers大小=" << feedComposition.size() << std::endl;
      assembleLocalJacobian(
        temperature, pressure, nV, vaporModel_, mV, muV);
      assembleLocalJacobian(
        temperature, pressure, nL, liquidModel_, mL, muL);
      
      std::cout<<"局部矩阵组装完成,迭代 "<<iter+1<<std::endl;

      // === NEW: 相内 Hessian 程序性校正（两相各做一次） ===
      double betaV = std::accumulate(nV.begin(), nV.end(), 0.0);
      double betaL = std::accumulate(nL.begin(), nL.end(), 0.0);
      std::vector<double> xV(nV.size()), xL(nL.size());
      for (size_t i=0;i<nV.size();++i){ xV[i] = (betaV>0)? nV[i]/betaV : 1.0/double(nV.size()); }
      for (size_t i=0;i<nL.size();++i){ xL[i] = (betaL>0)? nL[i]/betaL : 1.0/double(nL.size()); }

      PhaseFixOptions pfx;
      auto fixV = fix_phase_hessian_one_phase(xV, mV, pfx);
      if (fixV.applied){
        mV = fixV.m_fixed;
        std::cout << "[PhaseFix] Vapor: lam_min "<<fixV.lam_min_before
                  <<" -> "<<fixV.lam_min_after<<"\n";
      }
      auto fixL = fix_phase_hessian_one_phase(xL, mL, pfx);
      if (fixL.applied){
        mL = fixL.m_fixed;
        std::cout << "[PhaseFix] Liquid: lam_min "<<fixL.lam_min_before
                  <<" -> "<<fixL.lam_min_after<<"\n";
      }
      // 2) 反演 m_j -> M_j
      MV = invert(mV);
      ML = invert(mL);
      // for (size_t i = 0; i < C; ++i) {
      //   for (size_t k = 0; k < C; ++k) {
      //       std::cout << MV[i][k] << "  ";  // 同一行打印
      //   }
      //   std::cout << std::endl;  // 换行
      // }    

      std::cout<<"局部矩阵反演完成,迭代 "<<iter+1<<std::endl;

      // 3) 全局系统装配
      assembleGlobalSystem(
        temperature,
        MV, ML, muV, muL,
        nV, nL,
        elementMatrix,
        feedComposition, 
        Acoef, rhs
      );  
      std::cout<<"全局系统装配完成,迭代 "<<iter+1<<std::endl;
      // for(size_t i = 0; i < rhs.size(); ++i) {
      //   std::cout << "rhs[" << i << "] = " << rhs[i] << std::endl;
      // }

      // === 4) 解线性系统 ===
      double lin_resid = 0.0;
      sol = solveGlobalLinearSystem(Acoef, rhs, &lin_resid);
      // 这里保留原来的迭代号打印
      std::cout << "线性系统求解完成,迭代 " << iter+1 << std::endl;

      // 拆解解向量：[Δλ_1…Δλ_E, Δβ_v, Δβ_l]
      Lambda.assign(sol.begin(), sol.begin() + E);
      deltaBeta[0] = sol[E];
      deltaBeta[1] = sol[E+1];
      std::cout << "  Lambda: ";
      for (double dl : Lambda) std::cout << dl << " ";
      std::cout << "\n  deltaBeta_v=" << deltaBeta[0]
                << "  deltaBeta_l=" << deltaBeta[1] << std::endl;

      // === 5) 回代恢复 Δn 并输出 Δn ===
      // 保持原来的“更新前”块打印（你已有）
      {
        double betaV_dbg = std::accumulate(nV.begin(), nV.end(), 0.0);
        double betaL_dbg = std::accumulate(nL.begin(), nL.end(), 0.0);
        std::cout << "  更新前 : ";
        std::cout << "\n  nV: ";
        for (double v : nV) std::cout << v << " ";
        std::cout << "\n  nL: ";
        for (double v : nL) std::cout << v << " ";
        std::cout << "\n  betaV: " << betaV_dbg << ", betaL: " << betaL_dbg;
        std::cout << "\n  xV: ";
        for (double v : nV) std::cout << v/betaV_dbg << " ";
        std::cout << "\n  xL: ";
        for (double v : nL) std::cout << v/betaL_dbg << " ";
        std::cout << std::endl;
      }

      backSubstituteDeltas(
        temperature, elementMatrix, MV, ML, muV, muL,
        nV, nL, Lambda, deltaBeta, dnV, dnL);

      // === 6) 线搜索 + 更新 ===
      double alpha = applyUpdate(
        temperature, muV, muL, dnV, dnL, nV, nL);

      // === 7) 收敛判断 ===
      {
        auto conv = checkConvergence(
          iter + 1 , tol, muV, muL, elementMatrix, nV, nL, feedComposition);
        // 原逻辑：同时满足两个阈值才收敛
        if (conv.converged) {
          double betaV = std::accumulate(nV.begin(), nV.end(), 0.0);
          double betaL = std::accumulate(nL.begin(), nL.end(), 0.0);
          result.success           = true;
          result.vaporFraction     = betaV / (betaV + betaL);
          result.vaporComposition  = nV;
          result.liquidComposition = nL;
          result.iterations        = iter + 1;
          result.convergenceError  = conv.max_mu_diff;
          return result;
        }
      }
    }
  } catch (const std::exception& e){
    std::cerr << "[RandFlash] Exception caught during iterations: " << e.what() << std::endl;
    result.success = false;
    return result;
  }
  // 若未收敛
  result.success = false;
  return result;
}

