#include "linear_solver.hpp"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <stdexcept>

namespace ls {

namespace {

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

/// 从行主序 std::vector<double> 构造 Eigen::MatrixXd
inline Matrix makeEigenMatrix(int n, int m, const std::vector<double>& data)
{
  if (static_cast<int>(data.size()) != n * m) {
    throw std::invalid_argument("makeEigenMatrix: size mismatch");
  }
  Matrix M(n, m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      M(i, j) = data[i * m + j];
    }
  }
  return M;
}

/// 将 Eigen::MatrixXd 按行主序拷贝到 std::vector<double>
inline std::vector<double> toRowMajor(const Matrix& M)
{
  const int n = static_cast<int>(M.rows());
  const int m = static_cast<int>(M.cols());
  std::vector<double> out(n * m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      out[i * m + j] = M(i, j);
    }
  }
  return out;
}

} // unnamed namespace

class EigenSolver : public LinearSolverInterface {
public:
  std::vector<double> solveSPD(
      int n,
      const std::vector<double>& A,
      const std::vector<double>& b) override
  {
    if (static_cast<int>(b.size()) != n) {
      throw std::invalid_argument("solveSPD: RHS size mismatch");
    }
    if (static_cast<int>(A.size()) != n * n) {
      throw std::invalid_argument("solveSPD: matrix size mismatch");
    }

    Matrix matA = makeEigenMatrix(n, n, A);
    Vector vecB = Eigen::Map<const Vector>(b.data(), n);

    Eigen::LLT<Matrix> llt(matA);
    if (llt.info() != Eigen::Success) {
      throw std::runtime_error("solveSPD: Cholesky factorization failed");
    }
    Vector x = llt.solve(vecB);
    if (llt.info() != Eigen::Success) {
      throw std::runtime_error("solveSPD: solve failed");
    }

    std::vector<double> res(n);
    for (int i = 0; i < n; ++i) res[i] = x[i];
    return res;
  }

  std::vector<double> solveDense(
      int n,
      const std::vector<double>& A,
      const std::vector<double>& b,
      double* residual_norm = nullptr) override
  {
    if (static_cast<int>(b.size()) != n) {
      throw std::invalid_argument("solveDense: RHS size mismatch");
    }
    if (static_cast<int>(A.size()) != n * n) {
      throw std::invalid_argument("solveDense: matrix size mismatch");
    }

    Matrix matA = makeEigenMatrix(n, n, A);
    Vector vecB = Eigen::Map<const Vector>(b.data(), n);

    Eigen::JacobiSVD<Matrix> svd(
        matA, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Vector x = svd.solve(vecB);

    if (residual_norm) {
      Vector r = matA * x - vecB;
      *residual_norm = r.norm();
    }

    std::vector<double> res(n);
    for (int i = 0; i < n; ++i) res[i] = x[i];
    return res;
  }

  void eigenDecomposeSymmetric(
      int n,
      const std::vector<double>& A,
      std::vector<double>& eigenvalues,
      std::vector<double>& eigenvectors) override
  {
    if (static_cast<int>(A.size()) != n * n) {
      throw std::invalid_argument("eigenDecomposeSymmetric: matrix size mismatch");
    }

    Matrix matA = makeEigenMatrix(n, n, A);
    Eigen::SelfAdjointEigenSolver<Matrix> es(matA);
    if (es.info() != Eigen::Success) {
      throw std::runtime_error("eigenDecomposeSymmetric: eigen decomposition failed");
    }

    const Vector& evals = es.eigenvalues();
    const Matrix& evecs = es.eigenvectors(); // 列为特征向量

    eigenvalues.resize(n);
    for (int i = 0; i < n; ++i) {
      eigenvalues[i] = evals[i];
    }

    // 转为行主序、列为特征向量的存储方式
    eigenvectors.resize(n * n);
    for (int j = 0; j < n; ++j) {        // 列 j
      for (int i = 0; i < n; ++i) {      // 行 i
        eigenvectors[i * n + j] = evecs(i, j);
      }
    }
  }

  std::vector<double> invertSPD(
      int n,
      const std::vector<double>& A) override
  {
    if (static_cast<int>(A.size()) != n * n) {
      throw std::invalid_argument("invertSPD: matrix size mismatch");
    }

    Matrix matA = makeEigenMatrix(n, n, A);

    // Try Cholesky factorization first (fast for positive definite matrices)
    Eigen::LLT<Matrix> llt(matA);
    if (llt.info() == Eigen::Success) {
      Matrix I = Matrix::Identity(n, n);
      Matrix inv = llt.solve(I);
      if (llt.info() == Eigen::Success) {
        return toRowMajor(inv);
      }
    }

    // Fallback to SVD if Cholesky fails (more robust for near-singular matrices)
    Eigen::JacobiSVD<Matrix> svd(matA, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Matrix I = Matrix::Identity(n, n);
    Matrix inv = svd.solve(I);

    return toRowMajor(inv);
  }
};

std::unique_ptr<LinearSolverInterface> createEigenSolver()
{
  return std::make_unique<EigenSolver>();
}

} // namespace ls
