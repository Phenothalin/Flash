#pragma once
#include <vector>
#include <string>
#include "thermo_backend.hpp"
#include "phase_stability.hpp"
#include "RR_vl.hpp"

namespace rr_solver {

// RR求解选项
struct RRSolveOptions {
  bool enable_stability_analysis = false;  // 是否启用相稳定性分析
  double tolerance = 1e-6;                 // 收敛容差
  int max_outer_iterations = 50;          // 最大外层迭代次数
  ConvergenceMethod method = ConvergenceMethod::NEWTON_RAPHSON;  // 收敛方法
  phase_stability::StabilityOptions stability_opts;  // 稳定性分析选项
  bool verbose = false;                    // 是否输出详细信息
};

// RR闪蒸结果
struct RRFlashResult {
  // 相数
  int num_phases;  // 1, 2, 或 3

  // 两相结果
  double vapor_fraction;
  std::vector<double> vapor_composition;
  std::vector<double> liquid_composition;

  // 三相结果 (当num_phases == 3时有效)
  double vapor_fraction_vll;
  double liquid1_fraction;
  double liquid2_fraction;
  std::vector<double> vapor_composition_vll;
  std::vector<double> liquid1_composition;
  std::vector<double> liquid2_composition;

  // 收敛信息
  bool converged;
  int iterations;
  std::string message;

  RRFlashResult()
    : num_phases(0), vapor_fraction(0.0),
      vapor_fraction_vll(0.0), liquid1_fraction(0.0), liquid2_fraction(0.0),
      converged(false), iterations(0) {}
};

// 统一的PT闪蒸求解接口
RRFlashResult solvePTFlash(double pressure,
                           double temperature,
                           const std::vector<double>& z,
                           thermo::IThermoBackend& thermo,
                           const RRSolveOptions& options = RRSolveOptions());

// 统一的PV闪蒸求解接口
RRFlashResult solvePVFlash(double pressure,
                           double vapor_fraction,
                           const std::vector<double>& z,
                           thermo::IThermoBackend& thermo,
                           const RRSolveOptions& options = RRSolveOptions());

// 统一的TV闪蒸求解接口
RRFlashResult solveTVFlash(double temperature,
                           double vapor_fraction,
                           const std::vector<double>& z,
                           thermo::IThermoBackend& thermo,
                           const RRSolveOptions& options = RRSolveOptions());

} // namespace rr_solver
