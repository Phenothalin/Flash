#include "RR_vll.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace rr_vll {

VLLFlash::VLLFlash(double pressure,
                   double temperature,
                   const std::vector<double>& z,
                   thermo::IThermoBackend& thermo_backend,
                   const std::vector<double>& initial_vapor_comp,
                   const std::vector<double>& initial_liquid1_comp,
                   const std::vector<double>& initial_liquid2_comp)
  : P_(pressure), T_(temperature), z_(z), thermo_backend_(thermo_backend),
    beta_V_(0.0), beta_L1_(0.0) {

  size_t nc = z_.size();

  // 初始化组成向量
  if (!initial_vapor_comp.empty() && initial_vapor_comp.size() == nc) {
    y_ = initial_vapor_comp;
  } else {
    y_.resize(nc, 1.0 / nc);
  }

  if (!initial_liquid1_comp.empty() && initial_liquid1_comp.size() == nc) {
    x_L1_ = initial_liquid1_comp;
  } else {
    x_L1_.resize(nc, 1.0 / nc);
  }

  if (!initial_liquid2_comp.empty() && initial_liquid2_comp.size() == nc) {
    x_L2_ = initial_liquid2_comp;
  } else {
    x_L2_.resize(nc, 1.0 / nc);
  }

  K_V_.resize(nc, 1.0);
  K_L1_.resize(nc, 1.0);
}

void VLLFlash::initialize() {
  // 初始化相分率
  beta_V_ = 0.33;
  beta_L1_ = 0.33;

  size_t nc = z_.size();

  // 使用Wilson K值初始化气相K值
  thermo_backend_.wilsonK(T_, P_, K_V_);

  // 对于VLL体系（如水+烃类），两个液相应该是：
  // L1 (油相): 富含烃类
  // L2 (水相): 富含水
  //
  // K_L1,i = x_L1,i / x_L2,i
  // - 对于水: K_L1 << 1 (水主要在L2相)
  // - 对于烃类: K_L1 >> 1 (烃类主要在L1相)
  //
  // 识别策略：Wilson K值最小的组分通常是"水类"（极性、低挥发性）
  // 找到K值最小的组分作为"水类"组分

  size_t water_like_idx = 0;
  double min_K = K_V_[0];
  for (size_t i = 1; i < nc; ++i) {
    if (K_V_[i] < min_K) {
      min_K = K_V_[i];
      water_like_idx = i;
    }
  }

  // 初始化K_L1值，创建两个液相的明显分离
  for (size_t i = 0; i < nc; ++i) {
    if (i == water_like_idx) {
      // 水类组分：主要留在L2（水相），K_L1很小
      K_L1_[i] = 0.001;
    } else {
      // 烃类组分：主要进入L1（油相），K_L1很大
      // 使用K_V的倒数作为基础，确保烃类富集在L1
      K_L1_[i] = std::max(100.0, 1.0 / std::max(K_V_[i], 0.01));
    }
  }

  // 先计算一次相组成
  computePhaseCompositions();
}

void VLLFlash::computePhaseCompositions() {
  size_t nc = z_.size();
  double beta_L2 = 1.0 - beta_V_ - beta_L1_;

  // 从三相RR方程计算组成
  // D_i = 1 + β_V(K_V,i - 1) + β_L1(K_L1,i - 1)
  for (size_t i = 0; i < nc; ++i) {
    double D_i = 1.0 + beta_V_ * (K_V_[i] - 1.0) + beta_L1_ * (K_L1_[i] - 1.0);
    if (std::abs(D_i) < 1e-15) {
      throw std::runtime_error("VLL RR equation denominator near zero");
    }
    x_L2_[i] = z_[i] / D_i;
    y_[i] = K_V_[i] * x_L2_[i];
    x_L1_[i] = K_L1_[i] * x_L2_[i];
  }

  // 归一化
  double sum_y = 0.0, sum_L1 = 0.0, sum_L2 = 0.0;
  for (size_t i = 0; i < nc; ++i) {
    sum_y += y_[i];
    sum_L1 += x_L1_[i];
    sum_L2 += x_L2_[i];
  }

  for (size_t i = 0; i < nc; ++i) {
    y_[i] /= sum_y;
    x_L1_[i] /= sum_L1;
    x_L2_[i] /= sum_L2;
  }
}

void VLLFlash::updateKValues() {
  size_t nc = z_.size();

  // 计算各相的对数逸度系数
  auto lnphi_V = thermo_backend_.lnFugacityCoefficients(
      thermo::PhaseState{T_, P_, y_, thermo_backend_.vaporPhaseFlag()});

  auto lnphi_L1 = thermo_backend_.lnFugacityCoefficients(
      thermo::PhaseState{T_, P_, x_L1_, thermo_backend_.liquidPhaseFlag()});

  auto lnphi_L2 = thermo_backend_.lnFugacityCoefficients(
      thermo::PhaseState{T_, P_, x_L2_, thermo_backend_.liquidPhaseFlag()});

  // 更新K值: K_V,i = φ_L2,i / φ_V,i, K_L1,i = φ_L2,i / φ_L1,i
  // 使用对数形式避免数值问题: K = exp(ln(φ_L2) - ln(φ_other))
  for (size_t i = 0; i < nc; ++i) {
    K_V_[i] = std::exp(lnphi_L2[i] - lnphi_V[i]);
    K_L1_[i] = std::exp(lnphi_L2[i] - lnphi_L1[i]);
  }
}

bool VLLFlash::checkConvergence(const std::vector<double>& K_V_old,
                                const std::vector<double>& K_L1_old,
                                double beta_V_old,
                                double beta_L1_old,
                                double tolerance) {
  size_t nc = z_.size();

  // 检查K值变化
  double max_K_V_change = 0.0;
  double max_K_L1_change = 0.0;
  for (size_t i = 0; i < nc; ++i) {
    max_K_V_change = std::max(max_K_V_change, std::abs(K_V_[i] - K_V_old[i]));
    max_K_L1_change = std::max(max_K_L1_change, std::abs(K_L1_[i] - K_L1_old[i]));
  }

  // 检查相分率变化
  double beta_V_change = std::abs(beta_V_ - beta_V_old);
  double beta_L1_change = std::abs(beta_L1_ - beta_L1_old);

  return (max_K_V_change < tolerance && max_K_L1_change < tolerance &&
          beta_V_change < tolerance && beta_L1_change < tolerance);
}

bool VLLFlash::solveRRSystem(double tolerance, int max_iterations) {
  // 求解耦合的三相RR方程组:
  // f1(β_V, β_L1) = Σ z_i(K_V,i - 1) / D_i = 0
  // f2(β_V, β_L1) = Σ z_i(K_L1,i - 1) / D_i = 0
  // 其中 D_i = 1 + β_V(K_V,i - 1) + β_L1(K_L1,i - 1)

  size_t nc = z_.size();

  for (int iter = 0; iter < max_iterations; ++iter) {
    // 计算函数值
    double f1 = 0.0, f2 = 0.0;
    for (size_t i = 0; i < nc; ++i) {
      double D_i = 1.0 + beta_V_ * (K_V_[i] - 1.0) + beta_L1_ * (K_L1_[i] - 1.0);
      if (std::abs(D_i) < 1e-10) {
        D_i = D_i > 0 ? 1e-10 : -1e-10;  // 避免除零，但保持符号
      }
      f1 += z_[i] * (K_V_[i] - 1.0) / D_i;
      f2 += z_[i] * (K_L1_[i] - 1.0) / D_i;
    }

    // 检查收敛
    if (std::abs(f1) < tolerance && std::abs(f2) < tolerance) {
      return true;
    }

    // 计算Jacobian矩阵
    double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
    for (size_t i = 0; i < nc; ++i) {
      double K_V_m1 = K_V_[i] - 1.0;
      double K_L1_m1 = K_L1_[i] - 1.0;
      double D_i = 1.0 + beta_V_ * K_V_m1 + beta_L1_ * K_L1_m1;
      if (std::abs(D_i) < 1e-10) {
        D_i = D_i > 0 ? 1e-10 : -1e-10;
      }
      double D_i_sq = D_i * D_i;

      J11 += -z_[i] * K_V_m1 * K_V_m1 / D_i_sq;
      J12 += -z_[i] * K_V_m1 * K_L1_m1 / D_i_sq;
      J21 += -z_[i] * K_L1_m1 * K_V_m1 / D_i_sq;
      J22 += -z_[i] * K_L1_m1 * K_L1_m1 / D_i_sq;
    }

    // 求解线性系统: J * delta = -f
    double det = J11 * J22 - J12 * J21;

    double delta_beta_V, delta_beta_L1;

    if (std::abs(det) < 1e-12) {
      // Jacobian奇异，使用简单的梯度下降步
      double step = 0.1;
      delta_beta_V = -step * f1;
      delta_beta_L1 = -step * f2;
    } else {
      delta_beta_V = (-f1 * J22 + f2 * J12) / det;
      delta_beta_L1 = (f1 * J21 - f2 * J11) / det;
    }

    // 限制步长，避免过大的更新
    double max_step = 0.3;
    if (std::abs(delta_beta_V) > max_step) {
      delta_beta_V = delta_beta_V > 0 ? max_step : -max_step;
    }
    if (std::abs(delta_beta_L1) > max_step) {
      delta_beta_L1 = delta_beta_L1 > 0 ? max_step : -max_step;
    }

    // 更新相分率
    beta_V_ += delta_beta_V;
    beta_L1_ += delta_beta_L1;

    // 边界约束
    beta_V_ = std::max(0.001, std::min(0.998, beta_V_));
    beta_L1_ = std::max(0.001, std::min(0.998 - beta_V_, beta_L1_));
  }

  return false;
}

VLLResult VLLFlash::calculate(double tolerance, int max_iterations) {
  VLLResult result;

  // 初始化
  initialize();

  for (int outer_iter = 0; outer_iter < max_iterations; ++outer_iter) {
    // 保存旧值用于收敛检查
    std::vector<double> K_V_old = K_V_;
    std::vector<double> K_L1_old = K_L1_;
    double beta_V_old = beta_V_;
    double beta_L1_old = beta_L1_;

    // 1. 求解RR方程组得到新的相分率
    // 即使RR求解未完全收敛，也继续迭代
    solveRRSystem(tolerance, 50);

    // 2. 从相分率和K值计算相组成
    computePhaseCompositions();

    // 3. 从逸度系数更新K值
    updateKValues();

    // 4. 检查收敛
    if (checkConvergence(K_V_old, K_L1_old, beta_V_old, beta_L1_old, tolerance)) {
      result.converged = true;
      result.iterations = outer_iter + 1;
      result.vapor_fraction = beta_V_;
      result.liquid1_fraction = beta_L1_;
      result.liquid2_fraction = 1.0 - beta_V_ - beta_L1_;
      result.vapor_comp = y_;
      result.liquid1_comp = x_L1_;
      result.liquid2_comp = x_L2_;
      result.convergence_error = 0.0;
      return result;
    }
  }

  // 未收敛
  result.converged = false;
  result.iterations = max_iterations;
  result.vapor_fraction = beta_V_;
  result.liquid1_fraction = beta_L1_;
  result.liquid2_fraction = 1.0 - beta_V_ - beta_L1_;
  result.vapor_comp = y_;
  result.liquid1_comp = x_L1_;
  result.liquid2_comp = x_L2_;

  return result;
}

} // namespace rr_vll
