#include "RR_vll.hpp"
#include "thermo_backend.hpp"
#include "rr_solver.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <iostream>
#include <numeric>

// VLL三相闪蒸测试类
class VLLFlashTest : public ::testing::Test {
protected:
  // 辅助函数：验证物料平衡
  bool checkMassBalance(const std::vector<double>& z,
                        double beta_V, double beta_L1, double beta_L2,
                        const std::vector<double>& y,
                        const std::vector<double>& x_L1,
                        const std::vector<double>& x_L2,
                        double tolerance = 1e-4) {
    for (size_t i = 0; i < z.size(); ++i) {
      double calculated = beta_V * y[i] + beta_L1 * x_L1[i] + beta_L2 * x_L2[i];
      if (std::abs(calculated - z[i]) > tolerance) {
        std::cout << "Mass balance failed for component " << i
                  << ": z=" << z[i] << ", calculated=" << calculated
                  << ", diff=" << std::abs(calculated - z[i]) << "\n";
        return false;
      }
    }
    return true;
  }

  // 辅助函数：验证组成归一化
  bool checkNormalization(const std::vector<double>& composition,
                          double tolerance = 1e-6) {
    double sum = std::accumulate(composition.begin(), composition.end(), 0.0);
    return std::abs(sum - 1.0) < tolerance;
  }

  // 辅助函数：验证相分率之和为1
  bool checkPhaseFractionSum(double beta_V, double beta_L1, double beta_L2,
                             double tolerance = 1e-6) {
    double sum = beta_V + beta_L1 + beta_L2;
    return std::abs(sum - 1.0) < tolerance;
  }

  // 辅助函数：打印VLL结果
  void printVLLResult(const rr_vll::VLLResult& result,
                      const std::vector<std::string>& comp_names = {}) {
    std::cout << "\n--- VLL Flash Result ---\n";
    std::cout << "Converged: " << (result.converged ? "YES" : "NO") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Phase fractions:\n";
    std::cout << "  Vapor (V):   " << result.vapor_fraction << "\n";
    std::cout << "  Liquid1 (L1): " << result.liquid1_fraction << "\n";
    std::cout << "  Liquid2 (L2): " << result.liquid2_fraction << "\n";

    size_t nc = result.vapor_comp.size();
    std::cout << "\nPhase compositions:\n";
    std::cout << "Component\tVapor\t\tLiquid1\t\tLiquid2\n";
    for (size_t i = 0; i < nc; ++i) {
      if (i < comp_names.size()) {
        std::cout << comp_names[i] << "\t\t";
      } else {
        std::cout << "C" << i+1 << "\t\t";
      }
      std::cout << result.vapor_comp[i] << "\t"
                << result.liquid1_comp[i] << "\t"
                << result.liquid2_comp[i] << "\n";
    }
  }
};

// // 测试1: VLL基本功能测试 - 甲烷/乙烷/丙烷体系
// TEST_F(VLLFlashTest, VLL_Basic_C1C2C3) {
//   std::string components = "C1,C2,C3";
//   thermo::ThermoPackBackend thermo(components, "PR", "vdW", "Classic", "Default", false);

//   // 测试条件：低温高压，可能出现液液分相
//   double P = 5.0e6;  // 5 MPa
//   double T = 200.0;  // 200 K
//   std::vector<double> z = {0.3, 0.4, 0.3};

//   // 提供初始猜测
//   std::vector<double> init_vapor = {0.6, 0.3, 0.1};
//   std::vector<double> init_L1 = {0.2, 0.4, 0.4};
//   std::vector<double> init_L2 = {0.1, 0.3, 0.6};

//   rr_vll::VLLFlash flash(P, T, z, thermo, init_vapor, init_L1, init_L2);
//   auto result = flash.calculate(1e-6, 100);

//   std::cout << "\n=== VLL Basic Test (C1/C2/C3) ===\n";
//   std::cout << "T = " << T << " K, P = " << P/1e6 << " MPa\n";
//   printVLLResult(result, {"C1", "C2", "C3"});

//   // 验证基本约束
//   EXPECT_TRUE(checkNormalization(result.vapor_comp));
//   EXPECT_TRUE(checkNormalization(result.liquid1_comp));
//   EXPECT_TRUE(checkNormalization(result.liquid2_comp));
//   EXPECT_TRUE(checkPhaseFractionSum(result.vapor_fraction,
//                                      result.liquid1_fraction,
//                                      result.liquid2_fraction));

//   // 如果收敛，验证物料平衡
//   if (result.converged) {
//     EXPECT_TRUE(checkMassBalance(z, result.vapor_fraction,
//                                   result.liquid1_fraction,
//                                   result.liquid2_fraction,
//                                   result.vapor_comp,
//                                   result.liquid1_comp,
//                                   result.liquid2_comp));
//   }
// }

// // 测试2: 四组分体系VLL测试
// TEST_F(VLLFlashTest, VLL_FourComponent) {
//   std::string components = "C1,C2,C3,NC4";
//   thermo::ThermoPackBackend thermo(components, "PR");

//   double P = 4.0e6;  // 4 MPa
//   double T = 220.0;  // 220 K
//   std::vector<double> z = {0.25, 0.30, 0.25, 0.20};

//   rr_vll::VLLFlash flash(P, T, z, thermo);
//   auto result = flash.calculate(1e-6, 100);

//   std::cout << "\n=== VLL Four Component Test ===\n";
//   std::cout << "T = " << T << " K, P = " << P/1e6 << " MPa\n";
//   printVLLResult(result, {"C1", "C2", "C3", "nC4"});

//   EXPECT_TRUE(checkNormalization(result.vapor_comp));
//   EXPECT_TRUE(checkNormalization(result.liquid1_comp));
//   EXPECT_TRUE(checkNormalization(result.liquid2_comp));
//   EXPECT_TRUE(checkPhaseFractionSum(result.vapor_fraction,
//                                      result.liquid1_fraction,
//                                      result.liquid2_fraction));
// }

// // 测试3: 不同温度条件下的VLL测试
// TEST_F(VLLFlashTest, VLL_TemperatureVariation) {
//   std::string components = "C1,C2,C3";
//   thermo::ThermoPackBackend thermo(components, "PR");

//   double P = 3.0e6;  // 3 MPa
//   std::vector<double> z = {0.4, 0.35, 0.25};

//   // 测试多个温度点
//   std::vector<double> temperatures = {180.0, 200.0, 220.0};

//   for (double T : temperatures) {
//     std::cout << "\n=== VLL Temperature Test: T = " << T << " K ===\n";

//     rr_vll::VLLFlash flash(P, T, z, thermo);
//     auto result = flash.calculate(1e-6, 100);

//     std::cout << "Converged: " << (result.converged ? "YES" : "NO")
//               << ", Iterations: " << result.iterations << "\n";
//     std::cout << "Phase fractions: V=" << result.vapor_fraction
//               << ", L1=" << result.liquid1_fraction
//               << ", L2=" << result.liquid2_fraction << "\n";

//     EXPECT_TRUE(checkPhaseFractionSum(result.vapor_fraction,
//                                        result.liquid1_fraction,
//                                        result.liquid2_fraction));
//   }
// }

// // 测试4: 不同压力条件下的VLL测试
// TEST_F(VLLFlashTest, VLL_PressureVariation) {
//   std::string components = "C1,C2,C3";
//   thermo::ThermoPackBackend thermo(components, "PR");

//   double T = 200.0;  // 200 K
//   std::vector<double> z = {0.35, 0.35, 0.30};

//   std::vector<double> pressures = {2.0e6, 4.0e6, 6.0e6};

//   for (double P : pressures) {
//     std::cout << "\n=== VLL Pressure Test: P = " << P/1e6 << " MPa ===\n";

//     rr_vll::VLLFlash flash(P, T, z, thermo);
//     auto result = flash.calculate(1e-6, 100);

//     std::cout << "Converged: " << (result.converged ? "YES" : "NO")
//               << ", Iterations: " << result.iterations << "\n";
//     std::cout << "Phase fractions: V=" << result.vapor_fraction
//               << ", L1=" << result.liquid1_fraction
//               << ", L2=" << result.liquid2_fraction << "\n";

//     EXPECT_TRUE(checkPhaseFractionSum(result.vapor_fraction,
//                                        result.liquid1_fraction,
//                                        result.liquid2_fraction));
//   }
// }

// // 测试5: 统一求解器接口测试
// TEST_F(VLLFlashTest, VLL_UnifiedSolver) {
//   std::string components = "C1,C2,C3";
//   thermo::ThermoPackBackend thermo(components, "PR");

//   double P = 4.0e6;
//   double T = 210.0;
//   std::vector<double> z = {0.33, 0.34, 0.33};

//   rr_solver::RRSolveOptions options;
//   options.enable_stability_analysis = true;
//   options.verbose = true;

//   auto result = rr_solver::solvePTFlash(P, T, z, thermo, options);

//   std::cout << "\n=== Unified Solver Test ===\n";
//   std::cout << "Phases: " << result.num_phases << "\n";
//   std::cout << "Message: " << result.message << "\n";

//   EXPECT_TRUE(result.converged);
// }

// // 测试6: 五组分天然气体系
// TEST_F(VLLFlashTest, VLL_NaturalGas) {
//   std::string components = "C1,C2,C3,NC4,NC5";
//   thermo::ThermoPackBackend thermo(components, "PR");

//   double P = 5.0e6;
//   double T = 230.0;
//   std::vector<double> z = {0.50, 0.20, 0.15, 0.10, 0.05};

//   rr_vll::VLLFlash flash(P, T, z, thermo);
//   auto result = flash.calculate(1e-6, 100);

//   std::cout << "\n=== Natural Gas VLL Test ===\n";
//   printVLLResult(result, {"C1", "C2", "C3", "nC4", "nC5"});

//   EXPECT_TRUE(checkNormalization(result.vapor_comp));
//   EXPECT_TRUE(checkNormalization(result.liquid1_comp));
//   EXPECT_TRUE(checkNormalization(result.liquid2_comp));
// }

// // 测试7: 边界条件测试 - 极端组成
// TEST_F(VLLFlashTest, VLL_ExtremeComposition) {
//   std::string components = "C1,C2,C3";
//   thermo::ThermoPackBackend thermo(components, "PR");

//   double P = 3.0e6;
//   double T = 200.0;

//   // 轻组分为主
//   std::vector<double> z_light = {0.8, 0.15, 0.05};
//   rr_vll::VLLFlash flash_light(P, T, z_light, thermo);
//   auto result_light = flash_light.calculate(1e-6, 100);

//   std::cout << "\n=== Light-dominant composition ===\n";
//   std::cout << "V=" << result_light.vapor_fraction
//             << ", L1=" << result_light.liquid1_fraction
//             << ", L2=" << result_light.liquid2_fraction << "\n";

//   EXPECT_TRUE(checkPhaseFractionSum(result_light.vapor_fraction,
//                                      result_light.liquid1_fraction,
//                                      result_light.liquid2_fraction));
// }

// 测试8: 简单含水测试，与RNAD中的三相测试条件相同
TEST_F(VLLFlashTest, same_with_rand) {
  std::string components = "H2O,C1,nC6";
  thermo::ThermoPackBackend thermo(components, "SRK");

  double T = 298.15;     // Temperature in K
  double P = 1.01325e5;  // Pressure in Pa (1 atm)

  // 轻组分为主
  std::vector<double> z = {0.45, 0.05, 0.5};
  rr_vll::VLLFlash flash(P, T, z, thermo);
  auto result = flash.calculate(1e-6, 50);

  printVLLResult(result, {"H2O", "C1", "nC6"});  
  EXPECT_TRUE(checkPhaseFractionSum(result.vapor_fraction,
                                     result.liquid1_fraction,
                                     result.liquid2_fraction));
}
