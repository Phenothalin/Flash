#pragma once
#include <vector>
#include <memory>

namespace ls {

/// 对称正定线性系统求解接口： A x = b
struct LinearSolverInterface {
  virtual ~LinearSolverInterface() = default;
  virtual std::vector<double> solveSPD(
      int n,
      const std::vector<double>& A,
      const std::vector<double>& b) = 0;
};

/// 创建基于 Eigen 的解算器
std::unique_ptr<LinearSolverInterface> createEigenSolver();

} // namespace ls
