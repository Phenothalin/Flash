#pragma once
#include <vector>
#include <cmath>
#include <stdexcept>
#include <functional>

namespace rr_utils {

// RR方程工具类：提供Rachford-Rice方程的计算和求解方法
class RREquation {
public:
  // 计算RR方程值: f(β) = Σ z_i(K_i - 1) / [1 + β(K_i - 1)]
  static double evaluate(double beta,
                         const std::vector<double>& z,
                         const std::vector<double>& K);

  // 计算RR方程一阶导数: f'(β) = Σ -z_i(K_i - 1)^2 / [1 + β(K_i - 1)]^2
  static double derivative(double beta,
                           const std::vector<double>& z,
                           const std::vector<double>& K);

  // 计算RR方程二阶导数: f''(β) = Σ 2z_i(K_i - 1)^3 / [1 + β(K_i - 1)]^3
  static double secondDerivative(double beta,
                                 const std::vector<double>& z,
                                 const std::vector<double>& K);

  // 稳定形式的RR方程（当β > 0.5时使用）: f(β) = Σ z_i(1 - K_i) / [K_i + β(1 - K_i)]
  static double evaluateStable(double beta,
                               const std::vector<double>& z,
                               const std::vector<double>& K);

  // 从β和K值计算相组成
  static void computePhaseCompositions(double beta,
                                       const std::vector<double>& z,
                                       const std::vector<double>& K,
                                       std::vector<double>& x,  // 液相组成
                                       std::vector<double>& y); // 气相组成

  // Newton-Raphson求解器
  static double solveNewton(const std::vector<double>& z,
                           const std::vector<double>& K,
                           double beta_init,
                           double tolerance = 1e-10,
                           int max_iterations = 100);

  // Halley求解器（使用二阶导数）
  static double solveHalley(const std::vector<double>& z,
                           const std::vector<double>& K,
                           double beta_init,
                           double tolerance = 1e-10,
                           int max_iterations = 100);

  // 检查K值是否满足两相条件
  static bool checkTwoPhaseCondition(const std::vector<double>& z,
                                     const std::vector<double>& K);

private:
  // 辅助函数：计算分母 D_i = 1 + β(K_i - 1)
  static inline double denominator(double beta, double K_i) {
    return 1.0 + beta * (K_i - 1.0);
  }

  // 辅助函数：计算稳定形式的分母 D_i = K_i + β(1 - K_i)
  static inline double denominatorStable(double beta, double K_i) {
    return K_i + beta * (1.0 - K_i);
  }
};

} // namespace rr_utils
