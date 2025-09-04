#include "rand_flash.hpp"
#include "flash_calculator.hpp"
#include <Eigen/Dense>  
#include <numeric>
#include <cmath>
using namespace randflash;

static std::vector<std::vector<double>> invert(
  const std::vector<std::vector<double>>& mat)
{
  size_t n = mat.size();
  Eigen::MatrixXd M(n,n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      M(i,j) = mat[i][j];
  
  // 添加微小扰动避免奇异
  double jitter = 1e-10 * M.diagonal().mean();  // 基于对角线均值的扰动
  M.diagonal().array() += jitter;
  
  Eigen::MatrixXd Mi = M.inverse();
  std::vector<std::vector<double>> result(n, std::vector<double>(n));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      result[i][j] = Mi(i,j);  //将Eigen矩阵结果转换回二维vector格式
  return result;
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
  auto normalized = [](std::vector<double> v){
    double s = std::accumulate(v.begin(), v.end(), 0.0);
    if (s <= 0) throw std::invalid_argument("Composition sum must be positive");
    for (double &x : v) x /= s;
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
  // 1. 构造 PhaseState
  thermo::PhaseState state{ temperature, pressure, moleNumbers };

  // 2. 计算 μ 和 ∂μ/∂n
  mu = model.chemicalPotentials(state);
  for(size_t i = 0; i < moleNumbers.size(); ++i) {
    if (moleNumbers[i] <= 1e-10) {  // 允许微小正值,避免零
      moleNumbers[i] = 1e-10;  // 或抛出更友好的异常
      // throw std::invalid_argument("Mole number too small for component " + std::to_string(i));
    }
  }
  auto dmun = model.dMu_dN(state);
  // for(size_t i = 0; i < moleNumbers.size(); ++i) {
  //   for(size_t k = 0; k < moleNumbers.size(); ++k) {
  //     std::cout << "dmun[" << i << "][" << k << "] = "
  //               << dmun[i][k] << std::endl;
  //   }
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

      // 2) 反演 m_j -> M_j
      MV = invert(mV);
      ML = invert(mL);
    //   for(size_t i = 0; i < C; ++i) {
    //     for(size_t k = 0; k < C; ++k) {
    //       std::cout << "MV[" << i << "][" << k << "] = "
    //                 << MV[i][k] << std::endl;
    //     }
    //  }
      // for(size_t i = 0; i < C; ++i) {
      //   for(size_t k = 0; k < C; ++k) {
      //     std::cout << "ML[" << i << "][" << k << "] = "
      //               << ML[i][k] << std::endl;
      //   }
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

      // 4) 解线性系统 (E+2)x(E+2)
      int N = (int)rhs.size();
      // std::cout<<"Acoef size: "<<Acoef.size()<<", rhs size: "<<rhs.size()<<std::endl;
      assert((int)Acoef.size() == N*N);
      // Map the 1D Acoef into an Eigen matrix (row-major)
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
          A_mat(Acoef.data(), N, N);
          
      // 复制到一个本地变量,准备加 jitter
      Eigen::MatrixXd A_reg = A_mat;
      // **DEBUG**：打印原始全局矩阵
      // std::cout << "A_mat (size " << N << "×" << N << ") ***\n"
      //           << A_mat << "\n"
      //           << "*******************************\n";
      // 计算一个小的 regularization 强度：例如 trace(A) * 1e-8
      double traceA = A_reg.trace();
      double jitter = std::max(1e-8, traceA * 1e-6);
      // 在对角线上加 jitter
      A_reg.diagonal().array() += jitter;

      // std::cout << "*** DEBUG: A_reg (size " << N << "×" << N << ") ***\n"
      //           << A_reg << "\n"
      //           << "*******************************\n";

      // Map rhs into an Eigen vector
      Eigen::Map<Eigen::VectorXd> b_vec(rhs.data(), N);

      // Perform LLT (Cholesky) solve: A_mat * x = b_vec
      // Eigen::LLT<Eigen::MatrixXd> llt(A_reg);
      // if (llt.info() != Eigen::Success) {
      //   throw std::runtime_error("Cholesky decomposition failed even after regularization");
      // }

      // 用 SVD 求解最小二乘 / 最小范数解
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A_mat,
        Eigen::ComputeThinU | Eigen::ComputeThinV
      );
      Eigen::VectorXd x_vec = svd.solve(b_vec);

      // 检查残差
      double res_norm = (A_mat * x_vec - b_vec).norm();
      if (res_norm > 1e-6) {
      std::cerr << "[RandFlash] Warning: high residual = " << res_norm << "\n";
      }
    
      // 将结果写回 sol
      sol.assign(x_vec.data(), x_vec.data() + N);

      std::cout<<"线性系统求解完成,迭代 "<<iter+1<<std::endl;

      // 拆解解向量：[Δλ_1…Δλ_E, Δβ_v, Δβ_l]
      Lambda.assign(sol.begin(), sol.begin() + E);
      deltaBeta[0] = sol[E];
      deltaBeta[1] = sol[E+1];
      // 打印 Δλ 和 Δβ
      std::cout << "  Lambda: ";
      for (double dl : Lambda) std::cout << dl << " ";
      std::cout << "\n  deltaBeta_v=" << deltaBeta[0]
                << "  deltaBeta_l=" << deltaBeta[1] << std::endl;

      // 5) 回代恢复 Δn_{i,j}：
      double RT = R_CONST * temperature;
      double betaV = std::accumulate(nV.begin(), nV.end(), 0.0);
      double betaL = std::accumulate(nL.begin(), nL.end(), 0.0);

      // 在更新前,打印当前 nV, nL
      std::cout << "  更新前 : ";
      std::cout << "\n  nV: ";
      for (double v : nV) std::cout << v << " ";
      std::cout << "\n  nL: ";
      for (double v : nL) std::cout << v << " ";
      std::cout << "\n  betaV: " << betaV << ", betaL: " << betaL ;
      std::cout << "\n  xV: ";
      for (double v : nV) std::cout << v/betaV << " ";
      std::cout << "\n  xL: ";
      for (double v : nL) std::cout << v/betaL << " ";
      std::cout << std::endl;
      
      std::vector<double> lambdaComponent(C, 0.0);  // 组分维度的λ修正项
      for (size_t p = 0; p < C; ++p) {  // 遍历组分
        for (size_t ell = 0; ell < E; ++ell) {  // 遍历元素
          lambdaComponent[p] += elementMatrix[ell][p] * Lambda[ell];
        }
      }
      // 构造每相的 (A^T λ − μ/RT)
      std::vector<double> diffV(C), diffL(C);
      for (size_t p = 0; p < C; ++p) {
          diffV[p] = lambdaComponent[p] - muV[p] / RT;
          diffL[p] = lambdaComponent[p] - muL[p] / RT;
      }
      // 计算ΔnV和ΔnL
      for (size_t i = 0; i < C; ++i) {
        double xVi = nV[i] / betaV;   // 气相摩尔分数
        double xLi = nL[i] / betaL;   // 液相摩尔分数
    
        double combV = 0.0, combL = 0.0;
        for (size_t p = 0; p < C; ++p) {
            combV += MV[i][p] * diffV[p];
            combL += ML[i][p] * diffL[p];
        }
    
        dnV[i] = xVi * deltaBeta[0] + betaV * combV; 
        dnL[i] = xLi * deltaBeta[1] + betaL * combL; 
    }
      // 打印 Δn
      std::cout << "\n  dnV: ";
      for (double dv : dnV) std::cout << dv << " ";
      std::cout << "\n  dnL: ";
      for (double dv : dnL) std::cout << dv << " ";
      std::cout << "\n" <<std::endl;

      // 6) 线搜索,保持 n≥0
      std::vector<double> muV_over_RT(muV.size());
      std::vector<double> muL_over_RT(muL.size());
      for (size_t i = 0; i < muV.size(); ++i) muV_over_RT[i] = muV[i] / RT;
      for (size_t i = 0; i < muL.size(); ++i) muL_over_RT[i] = muL[i] / RT;

      // auto maxFeasibleAlpha = [&](const std::vector<double>& n, const std::vector<double>& dn){
      //   double a = 1.0;
      //   for(size_t i=0;i<n.size();++i){
      //       if (dn[i] < 0.0) a = std::min(a, 0.99 * n[i] / (-dn[i])); // 保证更新后仍为正
      //   }
      //   return a;
      // };
      // double alpha0 = std::min(1.0, std::min(maxFeasibleAlpha(nV, dnV),
      //                                      maxFeasibleAlpha(nL, dnL)));
      double alpha = lineSearch(nV, nL, dnV, dnL, muV_over_RT, muL_over_RT, 1, 1e-10, 0.5);
      std::cout << "  线搜索步长 alpha = " << alpha << std::endl;

      // 7) 更新
      for (size_t i = 0; i < C; ++i) {
        nV[i] += alpha * dnV[i];
        nL[i] += alpha * dnL[i];
      }
      betaV = std::accumulate(nV.begin(), nV.end(), 0.0);  // 气相总摩尔数 = 更新后各组分之和
      betaL = std::accumulate(nL.begin(), nL.end(), 0.0);  // 液相总摩尔数 = 更新后各组分之和
      
      // 打印更新后的数据（包含新的beta值）
      std::cout << "\n  更新后 : ";
      std::cout << "\n  nV: ";
      for (double v : nV) std::cout << v << " ";
      std::cout << "\n  nL: ";
      for (double v : nL) std::cout << v << " ";
      std::cout << "\n  betaV: " << betaV << ", betaL: " << betaL ;
      std::cout << "\n  xV: ";
      for (double v : nV) std::cout << (betaV > 1e-12 ? v/betaV : 0.0) << " ";
      std::cout << "\n  xL: ";
      for (double v : nL) std::cout << (betaL > 1e-12 ? v/betaL : 0.0) << " ";
      std::cout << "\n" <<std::endl;

      // 8) 收敛判据
      double max_mu_diff = 0.0;
      // 计算元素守恒残差（A×(nV+nL−feed) 的 L∞ 范数）
      auto elementResidualInf = [&](const std::vector<std::vector<double>>& A,
        const std::vector<double>& nV,
        const std::vector<double>& nL,
        const std::vector<double>& feed){
          const size_t E = A.size(), C = feed.size();
          double infnorm = 0.0;
          for (size_t ell=0; ell<E; ++ell){
            double acc = 0.0;
            for (size_t i=0; i<C; ++i){
              acc += A[ell][i] * ((nV[i]+nL[i]) - feed[i]);
            }
            infnorm = std::max(infnorm, std::abs(acc));
          }
          return infnorm;
      };
      double elem_err = elementResidualInf(elementMatrix, nV, nL, feedComposition);
      for (size_t i=0; i<C; ++i) {
        max_mu_diff = std::max(max_mu_diff, std::fabs(muV[i] - muL[i]));
      }
      std::cout << "Iter " << iter+1 << " | max_mu_diff: " << max_mu_diff << std::endl;
      std::cout << "Iter " << iter+1 << " | element error: " << elem_err << std::endl;
      if (max_mu_diff < tol && elem_err < 1e-8) {
        result.success           = true;
        result.vaporFraction     = betaV / (betaV + betaL);
        result.vaporComposition  = nV;
        result.liquidComposition = nL;
        result.iterations        = iter + 1;
        result.convergenceError  = max_mu_diff;
        return result;
      }
    }
  }catch (const std::exception& e) {
    std::cout << "捕获异常：" << e.what() << std::endl;
    return result; // 返回失败结果
  }

// 若未收敛
result.success = false;
return result;
}

