#include "rr_equation.hpp"
#include <algorithm>
#include <limits>

namespace rr_utils {

double RREquation::evaluate(double beta,
                            const std::vector<double>& z,
                            const std::vector<double>& K) {
  double sum = 0.0;
  for (size_t i = 0; i < z.size(); ++i) {
    double denom = denominator(beta, K[i]);
    if (std::abs(denom) < 1e-15) {
      throw std::runtime_error("RR equation denominator near zero");
    }
    sum += z[i] * (K[i] - 1.0) / denom;
  }
  return sum;
}

double RREquation::derivative(double beta,
                              const std::vector<double>& z,
                              const std::vector<double>& K) {
  double sum = 0.0;
  for (size_t i = 0; i < z.size(); ++i) {
    double Ki_minus_1 = K[i] - 1.0;
    double denom = denominator(beta, K[i]);
    if (std::abs(denom) < 1e-15) {
      throw std::runtime_error("RR equation denominator near zero");
    }
    sum += -z[i] * Ki_minus_1 * Ki_minus_1 / (denom * denom);
  }
  return sum;
}

double RREquation::secondDerivative(double beta,
                                    const std::vector<double>& z,
                                    const std::vector<double>& K) {
  double sum = 0.0;
  for (size_t i = 0; i < z.size(); ++i) {
    double Ki_minus_1 = K[i] - 1.0;
    double denom = denominator(beta, K[i]);
    if (std::abs(denom) < 1e-15) {
      throw std::runtime_error("RR equation denominator near zero");
    }
    sum += 2.0 * z[i] * Ki_minus_1 * Ki_minus_1 * Ki_minus_1 /
           (denom * denom * denom);
  }
  return sum;
}

double RREquation::evaluateStable(double beta,
                                  const std::vector<double>& z,
                                  const std::vector<double>& K) {
  double sum = 0.0;
  for (size_t i = 0; i < z.size(); ++i) {
    double denom = denominatorStable(beta, K[i]);
    if (std::abs(denom) < 1e-15) {
      throw std::runtime_error("RR equation denominator near zero");
    }
    sum += z[i] * (1.0 - K[i]) / denom;
  }
  return sum;
}

void RREquation::computePhaseCompositions(double beta,
                                          const std::vector<double>& z,
                                          const std::vector<double>& K,
                                          std::vector<double>& x,
                                          std::vector<double>& y) {
  size_t nc = z.size();
  x.resize(nc);
  y.resize(nc);

  double sum_x = 0.0;
  double sum_y = 0.0;

  for (size_t i = 0; i < nc; ++i) {
    double denom = denominator(beta, K[i]);
    if (std::abs(denom) < 1e-15) {
      throw std::runtime_error("RR equation denominator near zero");
    }
    x[i] = z[i] / denom;
    y[i] = K[i] * x[i];
    sum_x += x[i];
    sum_y += y[i];
  }

  // 归一化
  for (size_t i = 0; i < nc; ++i) {
    x[i] /= sum_x;
    y[i] /= sum_y;
  }
}

bool RREquation::checkTwoPhaseCondition(const std::vector<double>& z,
                                        const std::vector<double>& K) {
  // 检查是否满足两相条件：min(K) < 1 < max(K)
  double K_min = *std::min_element(K.begin(), K.end());
  double K_max = *std::max_element(K.begin(), K.end());

  return (K_min < 1.0 && K_max > 1.0);
}

double RREquation::solveNewton(const std::vector<double>& z,
                               const std::vector<double>& K,
                               double beta_init,
                               double tolerance,
                               int max_iterations) {
  double beta = beta_init;

  // 确保初始值在有效范围内
  beta = std::max(0.0, std::min(1.0, beta));

  for (int iter = 0; iter < max_iterations; ++iter) {
    double f = evaluate(beta, z, K);
    double df = derivative(beta, z, K);

    if (std::abs(df) < 1e-15) {
      // 导数接近零，尝试使用稳定形式
      if (beta > 0.5) {
        f = evaluateStable(beta, z, K);
      }
      if (std::abs(df) < 1e-15) {
        break; // 无法继续迭代
      }
    }

    double delta = -f / df;
    double beta_new = beta + delta;

    // 边界约束
    beta_new = std::max(0.0, std::min(1.0, beta_new));

    if (std::abs(beta_new - beta) < tolerance) {
      return beta_new;
    }

    beta = beta_new;
  }

  return beta;
}

double RREquation::solveHalley(const std::vector<double>& z,
                               const std::vector<double>& K,
                               double beta_init,
                               double tolerance,
                               int max_iterations) {
  double beta = beta_init;

  // 确保初始值在有效范围内
  beta = std::max(0.0, std::min(1.0, beta));

  for (int iter = 0; iter < max_iterations; ++iter) {
    double f = evaluate(beta, z, K);
    double df = derivative(beta, z, K);
    double d2f = secondDerivative(beta, z, K);

    if (std::abs(df) < 1e-15) {
      break; // 无法继续迭代
    }

    // Halley's method: x_new = x - 2f*f' / (2(f')^2 - f*f'')
    double numerator = 2.0 * f * df;
    double denominator = 2.0 * df * df - f * d2f;

    if (std::abs(denominator) < 1e-15) {
      // 退化为Newton法
      double delta = -f / df;
      beta = beta + delta;
    } else {
      double delta = -numerator / denominator;
      beta = beta + delta;
    }

    // 边界约束
    beta = std::max(0.0, std::min(1.0, beta));

    if (std::abs(f) < tolerance) {
      return beta;
    }
  }

  return beta;
}

} // namespace rr_utils
