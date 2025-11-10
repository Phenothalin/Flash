#pragma once
#include <vector>
#include <memory>
#include <Eigen/Dense>  

namespace ls {

class EigenSolver : public LinearSolverInterface {
  public:
    std::vector<double> solveSPD(
        int n,
        const std::vector<double>& A,
        const std::vector<double>& b) override
    {
      if ((int)b.size() != n) throw std::invalid_argument("RHS size mismatch");
      Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>
        matA(A.data(), n, n);
      Eigen::Map<const Eigen::VectorXd> vecB(b.data(), n);
  
      Eigen::LLT<Eigen::MatrixXd> llt(matA);
      if (llt.info() != Eigen::Success) throw std::runtime_error("Cholesky failed");
      Eigen::VectorXd sol = llt.solve(vecB);
      if (llt.info() != Eigen::Success) throw std::runtime_error("Solve failed");
  
      return std::vector<double>(sol.data(), sol.data() + n);
    }
  };
  
std::unique_ptr<LinearSolverInterface> createEigenSolver() {
  return std::make_unique<EigenSolver>();
} 

/// 对称正定线性系统求解接口： A x = b
struct LinearSolverInterface {
  virtual ~LinearSolverInterface() = default;
  virtual std::vector<double> solveSPD(
      int n,
      const std::vector<double>& A,
      const std::vector<double>& b) = 0;
};


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

// === utils: std::vector <-> Eigen ===
static Eigen::VectorXd toEig(const std::vector<double>& v){
  Eigen::VectorXd e(v.size());
  for (size_t i=0;i<v.size();++i) e[i]=v[i];
  return e;
}
static Eigen::MatrixXd toEig(const std::vector<std::vector<double>>& M){
  const size_t r=M.size(), c=M[0].size();
  Eigen::MatrixXd E(r,c);
  for (size_t i=0;i<r;++i) for (size_t j=0;j<c;++j) E(i,j)=M[i][j];
  return E;
}
static std::vector<std::vector<double>> toStd(const Eigen::MatrixXd& E){
  std::vector<std::vector<double>> M(E.rows(), std::vector<double>(E.cols()));
  for (int i=0;i<E.rows();++i) for (int j=0;j<E.cols();++j) M[i][j]=E(i,j);
  return M;
}
} // namespace ls
