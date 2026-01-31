#include "rr_solver.hpp"
#include "RR_vll.hpp"
#include <iostream>

namespace rr_solver {

RRFlashResult solvePTFlash(double pressure,
                           double temperature,
                           const std::vector<double>& z,
                           thermo::IThermoBackend& thermo,
                           const RRSolveOptions& options) {
  RRFlashResult result;

  if (options.verbose) {
    std::cout << "\n=== RR PT Flash Solver ===\n";
    std::cout << "P = " << pressure/1e6 << " MPa, T = " << temperature << " K\n";
    std::cout << "Stability analysis: " << (options.enable_stability_analysis ? "ON" : "OFF") << "\n";
  }

  // 如果启用稳定性分析
  if (options.enable_stability_analysis) {
    phase_stability::PhaseStabilityAnalyzer analyzer(thermo);
    auto stability_result = analyzer.analyze(temperature, pressure, z, options.stability_opts);

    if (options.verbose) {
      std::cout << "Stability analysis result: "
                << (stability_result.stable ? "STABLE" : "UNSTABLE") << "\n";
      std::cout << "Number of incipient phases: " << stability_result.incipient.size() << "\n";
    }

    // 单相稳定
    if (stability_result.stable) {
      result.num_phases = 1;
      result.converged = true;
      result.vapor_fraction = 1.0;  // 或根据相态判断
      result.vapor_composition = z;
      result.message = "Single phase stable";
      return result;
    }

    // 检测到多个新相 - 可能是三相
    if (stability_result.incipient.size() >= 2) {
      if (options.verbose) {
        std::cout << "Attempting three-phase VLL flash...\n";
      }

      // 使用检测到的相组成作为初始猜测
      std::vector<double> init_vapor, init_liquid1, init_liquid2;

      // 简单策略：使用前两个新相作为初始猜测
      if (stability_result.incipient.size() >= 2) {
        init_liquid1 = stability_result.incipient[0].x;
        init_liquid2 = stability_result.incipient[1].x;
      }

      try {
        rr_vll::VLLFlash vll_flash(pressure, temperature, z, thermo,
                                    init_vapor, init_liquid1, init_liquid2);
        auto vll_result = vll_flash.calculate(options.tolerance, options.max_outer_iterations);

        if (vll_result.converged) {
          result.num_phases = 3;
          result.converged = true;
          result.vapor_fraction_vll = vll_result.vapor_fraction;
          result.liquid1_fraction = vll_result.liquid1_fraction;
          result.liquid2_fraction = vll_result.liquid2_fraction;
          result.vapor_composition_vll = vll_result.vapor_comp;
          result.liquid1_composition = vll_result.liquid1_comp;
          result.liquid2_composition = vll_result.liquid2_comp;
          result.iterations = vll_result.iterations;
          result.message = "Three-phase VLL equilibrium";
          return result;
        }
      } catch (const std::exception& e) {
        if (options.verbose) {
          std::cout << "VLL flash failed: " << e.what() << "\n";
          std::cout << "Falling back to two-phase flash...\n";
        }
      }
    }
  }

  // 执行标准两相闪蒸
  if (options.verbose) {
    std::cout << "Performing two-phase flash...\n";
  }

  PTFlash flash(pressure, temperature, z, thermo);
  flash.calculate(options.method);

  result.num_phases = 2;
  result.converged = true;
  result.vapor_fraction = flash.getVaporFraction();
  result.vapor_composition = flash.getVapComp();
  result.liquid_composition = flash.getLiqComp();
  result.message = "Two-phase equilibrium";

  return result;
}

RRFlashResult solvePVFlash(double pressure,
                           double vapor_fraction,
                           const std::vector<double>& z,
                           thermo::IThermoBackend& thermo,
                           const RRSolveOptions& options) {
  RRFlashResult result;

  if (options.verbose) {
    std::cout << "\n=== RR PV Flash Solver ===\n";
    std::cout << "P = " << pressure/1e6 << " MPa, V = " << vapor_fraction << "\n";
  }

  // PV闪蒸通常不需要稳定性分析（已知相分率）
  PVFlash flash(pressure, vapor_fraction, z, thermo);
  flash.calculate(options.method);

  result.num_phases = 2;
  result.converged = true;
  result.vapor_fraction = flash.getVaporFraction();
  result.vapor_composition = flash.getVapComp();
  result.liquid_composition = flash.getLiqComp();
  result.message = "Two-phase PV flash";

  return result;
}

RRFlashResult solveTVFlash(double temperature,
                           double vapor_fraction,
                           const std::vector<double>& z,
                           thermo::IThermoBackend& thermo,
                           const RRSolveOptions& options) {
  RRFlashResult result;

  if (options.verbose) {
    std::cout << "\n=== RR TV Flash Solver ===\n";
    std::cout << "T = " << temperature << " K, V = " << vapor_fraction << "\n";
  }

  // TV闪蒸通常不需要稳定性分析（已知相分率）
  TVFlash flash(temperature, vapor_fraction, z, thermo);
  flash.calculate(options.method);

  result.num_phases = 2;
  result.converged = true;
  result.vapor_fraction = flash.getVaporFraction();
  result.vapor_composition = flash.getVapComp();
  result.liquid_composition = flash.getLiqComp();
  result.message = "Two-phase TV flash";

  return result;
}

} // namespace rr_solver
