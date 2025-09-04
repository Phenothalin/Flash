#pragma once
#include <vector>
#include <memory>
#include "property_package.hpp"

static constexpr double R_CONST = 8.314462618;

namespace thermo {

/// 描述一个相态的温度、压力和各组分摩尔量
struct PhaseState {
  double temperature;                ///< 温度 (K)
  double pressure;                   ///< 压力 (Pa)
  std::vector<double> moleNumbers;   ///< 各组分的绝对摩尔数
};

/// 化学势模型的抽象接口
class ThermoModel {
public:
  virtual ~ThermoModel() = default;

  /// 给定相态，返回每个组分的化学势 μ_i
  virtual std::vector<double> chemicalPotentials(const PhaseState& state) = 0;

  /// 返回 ∂μ_i/∂n_k（Hessian 块矩阵）
  virtual std::vector<std::vector<double>> dMu_dN(const PhaseState& state) = 0;

  /// 计算 dPhi_dn 矩阵（逸度系数对摩尔数的导数）
  virtual std::vector<std::vector<double>> dPhi_dn(const PhaseState& state) = 0;
  /// 温度、压力导数（可选）
  virtual std::vector<double> dMu_dT(const PhaseState& state) = 0;
  virtual std::vector<double> dMu_dP(const PhaseState& state) = 0;
};

/// 现有 property_package 适配到 ThermoModel
class PropertyPackageAdapter : public ThermoModel {
public:
  /// 构造：ppType 指定状态方程类型，componentCluster 指定组分信息
  PropertyPackageAdapter(
    PropertyPackageType ppType,
    const std::shared_ptr<material_object::Cluster>& componentCluster)
    : prop_(ppType, componentCluster){}

  /// 设置当前求解的相：0 = vapor, 1 = liquid
  void setPhaseIndex(int index) { phaseIndex_ = index; }

  
  /// 化学势 μ_i = R T (ln φ_i + ln x_i)
  std::vector<double> chemicalPotentials(const PhaseState& state) override {
    size_t componentCount = state.moleNumbers.size();
    // 1) 归一化出摩尔分率
    double total = 0.0;
    for (auto n : state.moleNumbers) total += n;
    std::vector<double> molefraction(componentCount);
    for (size_t i = 0; i < componentCount; ++i) molefraction[i] = state.moleNumbers[i] / total;

    auto fugacityCoefficients = (phaseIndex_ == 0)
      ? prop_.calculateVaporFugacityCoefficientMixture(
          state.temperature, state.pressure, molefraction)
      : prop_.calculateLiquidFugacityCoefficientMixture(
          state.temperature, state.pressure, molefraction);
    std::vector<double> chemicalPotent(componentCount);
    for (size_t i = 0; i < componentCount; ++i) {
      chemicalPotent[i] = R_CONST * state.temperature * (
        std::log(fugacityCoefficients[i]) +
        std::log(molefraction[i]) + std::log(state.pressure)
      );
      // std::cout << "Fugacity coefficient: "
      // << fugacityCoefficients[i] << std::endl;
      // std::cout << "Component " << i << ": "
      //           << "μ = " << chemicalPotent[i] << " J/mol" << std::endl;
    }
    return chemicalPotent;
  }
  /// 公开的 dPhi_dn 计算接口（外部可直接调用）
  std::vector<std::vector<double>> dPhi_dn(const PhaseState& state) override {
    size_t componentCount = state.moleNumbers.size();
    // 1) 计算总摩尔数和摩尔分率
    double total = 0.0;
    for (auto n : state.moleNumbers) total += n;
    std::vector<double> molefraction(componentCount);
    for (size_t i = 0; i < componentCount; ++i) {
      molefraction[i] = state.moleNumbers[i] / total;
    }

    // 2) 获取逸度系数对组成的导数 dPhi_dX
    auto dPhi_dX = (phaseIndex_ == 0)
      ? prop_.calculateVaporFugacityCoefficientMixtureDerivativeComposition(
          state.temperature, state.pressure, molefraction)
      : prop_.calculateLiquidFugacityCoefficientMixtureDerivativeComposition(
          state.temperature, state.pressure, molefraction); 

    // 3) 计算 dPhi_dn 矩阵
    std::vector<std::vector<double>> dPhi_dn_matrix(componentCount,
      std::vector<double>(componentCount, 0.0));
    
    for (size_t i = 0; i < componentCount; ++i) {
      for (size_t k = 0; k < componentCount; ++k) {
        for (size_t l = 0; l < componentCount; ++l) {
          double dXl_dn = (l == k
            ? (total - state.moleNumbers[l])
            : -state.moleNumbers[l]
          ) / (total * total);
          dPhi_dn_matrix[i][k] += dPhi_dX[i][l] * dXl_dn;
        }
      }
    }

    return dPhi_dn_matrix;
  }
  /// ∂μ_i/∂n_k
  std::vector<std::vector<double>> dMu_dN(const PhaseState& state) override {
    size_t componentCount = state.moleNumbers.size();
    // 1) 归一化出摩尔分率
    double total = 0.0;
    for (auto n : state.moleNumbers) total += n;
    std::vector<double> molefraction(componentCount);
    for (size_t i = 0; i < componentCount; ++i) molefraction[i] = state.moleNumbers[i] / total;
    
    auto fugacityCoefficients = (phaseIndex_ == 0)
      ? prop_.calculateVaporFugacityCoefficientMixture(
          state.temperature, state.pressure, molefraction)
      : prop_.calculateLiquidFugacityCoefficientMixture(
          state.temperature, state.pressure, molefraction);
    auto dPhi_dX = (phaseIndex_ == 0)
      ? prop_.calculateVaporFugacityCoefficientMixtureDerivativeComposition(
          state.temperature, state.pressure, molefraction)
      : prop_.calculateLiquidFugacityCoefficientMixtureDerivativeComposition(
          state.temperature, state.pressure, molefraction); 

    // 计算 dPhi_dn 矩阵（调用独立函数）
    auto dPhi_dn_matrix = dPhi_dn(state);

    // 打印 dPhi_dn 矩阵
    // std::cout << "dPhi_dn matrix:" << std::endl;
    // for (size_t i = 0; i < componentCount; ++i) {
    //   for (size_t k = 0; k < componentCount; ++k) {
    //     std::cout << "dPhi_dn[" << i << "][" << k << "] = "
    //               << dPhi_dn[i][k] << std::endl;
    //   }
    // }

    // 计算 Jacobian 矩阵
    std::vector<std::vector<double>> jacobian(componentCount,
      std::vector<double>(componentCount, 0.0));
    for (size_t i = 0; i < componentCount; ++i) {
      for (size_t k = 0; k < componentCount; ++k) {
        double term = (dPhi_dn_matrix[i][k] / fugacityCoefficients[i])
          + (i == k ? 1.0 / state.moleNumbers[i] : 0.0)
          - 1.0 / total;
        jacobian[i][k] = R_CONST * state.temperature * term;
      }
    }
    return jacobian;
  }

  std::vector<double> dMu_dT(const PhaseState&) override {
    return {};
  }
  std::vector<double> dMu_dP(const PhaseState&) override {
    return {};
  }

private:
  property_package::PropertyPackage prop_;
  int phaseIndex_ = 0;
};

} // namespace thermo
