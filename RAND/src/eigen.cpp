#include "linear_solver.hpp"
#include <Eigen/Dense>
#include <stdexcept>

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

} // namespace ls
