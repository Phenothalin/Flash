#include "RR_vl.hpp"
#include "thermo_backend.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <iostream>

// 基础测试类
class RRFlashTest : public ::testing::Test {
protected:
  void SetUp() override {
    // 测试前准备
  }

  void TearDown() override {
    // 测试后清理
  }

  // 辅助函数：验证物料平衡
  bool checkMassBalance(const std::vector<double>& z,
                        double beta,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        double tolerance = 1e-5) {
    for (size_t i = 0; i < z.size(); ++i) {
      double calculated = (1.0 - beta) * x[i] + beta * y[i];
      if (std::abs(calculated - z[i]) > tolerance) {
        std::cout << "Mass balance failed for component " << i
                  << ": z=" << z[i] << ", calculated=" << calculated
                  << ", error=" << std::abs(calculated - z[i]) << "\n";
        return false;
      }
    }
    return true;
  }

  // 辅助函数：验证组成归一化
  bool checkNormalization(const std::vector<double>& composition,
                          double tolerance = 1e-6) {
    double sum = 0.0;
    for (double x : composition) {
      sum += x;
    }
    return std::abs(sum - 1.0) < tolerance;
  }
};

// // 测试1: 基本两相PT闪蒸 - 甲烷/乙烷体系
// TEST_F(RRFlashTest, PTFlash_Basic) {
//   // 创建 ThermoPackBackend (PR方程)
//   std::string components = "C1,C2";  // 甲烷, 乙烷
//   thermo::ThermoPackBackend thermo(components, "PR", "vdW", "Classic", "Default", false);

//   // 测试条件
//   double P = 2.0e6;  // 2 MPa
//   double T = 250.0;  // 250 K
//   std::vector<double> z = {0.5, 0.5};  // 等摩尔混合物

//   // 创建PTFlash对象
//   PTFlash flash(P, T, z, thermo);

//   // 执行计算
//   flash.calculate(ConvergenceMethod::NEWTON_RAPHSON);

//   // 获取结果
//   double beta = flash.getVaporFraction();
//   std::vector<double> y = flash.getVapComp();
//   std::vector<double> x = flash.getLiqComp();

//   // 输出结果
//   std::cout << "\n=== PTFlash Basic Test ===\n";
//   std::cout << "T = " << T << " K, P = " << P/1e6 << " MPa\n";
//   std::cout << "Vapor fraction: " << beta << "\n";
//   std::cout << "Vapor composition: [" << y[0] << ", " << y[1] << "]\n";
//   std::cout << "Liquid composition: [" << x[0] << ", " << x[1] << "]\n";

//   // 验证结果
//   EXPECT_GE(beta, 0.0);
//   EXPECT_LE(beta, 1.0);
//   EXPECT_TRUE(checkNormalization(y));
//   EXPECT_TRUE(checkNormalization(x));
//   EXPECT_TRUE(checkMassBalance(z, beta, x, y));

//   // 验证轻组分在气相中富集
//   EXPECT_GT(y[0], x[0]);  // 甲烷在气相中更多
// }

// // 测试2: 三组分体系PT闪蒸
// TEST_F(RRFlashTest, PTFlash_ThreeComponent) {
//   // 甲烷/乙烷/丙烷体系
//   std::string components = "C1,C2,C3";
//   thermo::ThermoPackBackend thermo(components, "PR", "vdW", "Classic", "Default", false);

//   double P = 3.0e6;  // 3 MPa
//   double T = 280.0;  // 280 K
//   std::vector<double> z = {0.4, 0.35, 0.25};

//   PTFlash flash(P, T, z, thermo);
//   flash.calculate(ConvergenceMethod::NEWTON_RAPHSON);

//   double beta = flash.getVaporFraction();
//   std::vector<double> y = flash.getVapComp();
//   std::vector<double> x = flash.getLiqComp();

//   std::cout << "\n=== PTFlash Three Component Test ===\n";
//   std::cout << "T = " << T << " K, P = " << P/1e6 << " MPa\n";
//   std::cout << "Vapor fraction: " << beta << "\n";
//   std::cout << "Vapor: [" << y[0] << ", " << y[1] << ", " << y[2] << "]\n";
//   std::cout << "Liquid: [" << x[0] << ", " << x[1] << ", " << x[2] << "]\n";

//   EXPECT_GE(beta, 0.0);
//   EXPECT_LE(beta, 1.0);
//   EXPECT_TRUE(checkNormalization(y));
//   EXPECT_TRUE(checkNormalization(x));
//   EXPECT_TRUE(checkMassBalance(z, beta, x, y));
// }

// 测试3: 复杂多组分体系PT闪蒸
TEST_F(RRFlashTest, PTFlash_MultiComponent) {
      std::string comps =
      "N2,CO2,C1,C2,C3,"
      "iC4,nC4,iC5,nC5,nC6,nC7";

    thermo::ThermoPackBackend backend(
      comps, "SRK", "vdW", "Classic", "Default", false);

    double T = 295.0;      // K
    double P = 2.0e6;      // Pa
    std::vector<double> feed = {
      0.003, 0.015, 0.55, 0.14, 0.12,
      0.05, 0.045, 0.03, 0.025, 0.012, 0.01
    };

    PTFlash flash(P, T, feed, backend);
    flash.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    double beta = flash.getVaporFraction();
    std::vector<double> y = flash.getVapComp();
    std::vector<double> x = flash.getLiqComp();
    std::cout << "\n=== PTFlash Multi-Component Test ===\n";
    std::cout << "T = " << T << " K, P = " << P/1e6 << " MPa\n";
    std::cout << "Vapor fraction: " << beta << "\n";
    std::cout << "Vapor composition: ";
    for (double comp : y) {
      std::cout << comp << " ";
    }
    std::cout << "\nLiquid composition: ";
    for (double comp : x) {
      std::cout << comp << " ";
    }
    std::cout << "\n";  
  }
