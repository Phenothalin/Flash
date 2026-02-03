// RAND/tests/stable_method_test.cpp
// 测试稳定法相数判断功能

#include <gtest/gtest.h>
#include "rand_flash.hpp"
#include "thermopack_adapter.hpp"
#include "linear_solver.hpp"
#include <vector>
#include <cmath>
#include <iostream>

using namespace randflash;
using namespace thermo;

class StableMethodTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 创建热力学后端和线性求解器
        backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2,C3", "PR");
        linSolver = ls::createEigenSolver();
        flash = std::make_unique<RandFlash>(*backend, *linSolver);

        // 单位元素矩阵（非反应体系）
        elementMatrix = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };
    }

    std::unique_ptr<thermo::IThermoBackend> backend;
    std::shared_ptr<ls::LinearSolverInterface> linSolver;
    std::unique_ptr<RandFlash> flash;
    std::vector<std::vector<double>> elementMatrix;
};

// Test 1: 单相稳定系统（超临界流体）
TEST_F(StableMethodTest, SinglePhaseStable) {
    FlashInput input;
    input.temperature = 400.0;  // K - 高温
    input.pressure = 1e7;       // Pa - 高压（超临界）
    input.feedMoles = {0.4, 0.3, 0.3};

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 5;
    opt.max_phases = 3;
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, elementMatrix, opt);

    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.numPhases(), 1) << "超临界条件应为单相";

    std::cout << "Test 1 - Single Phase Stable: "
              << result.numPhases() << " phases, "
              << "iterations: " << result.iterations << std::endl;
}

// Test 2: 简单VLE（二元混合物）
TEST_F(StableMethodTest, SimpleVLE) {
    FlashInput input;
    input.temperature = 250.0;  // K - 中等温度
    input.pressure = 2e6;       // Pa - 中等压力
    input.feedMoles = {0.5, 0.3, 0.2};

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 10;
    opt.max_phases = 3;
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, elementMatrix, opt);

    EXPECT_TRUE(result.success);
    EXPECT_GE(result.numPhases(), 1);
    EXPECT_LE(result.numPhases(), 3);

    std::cout << "Test 2 - Simple VLE: "
              << result.numPhases() << " phases, "
              << "iterations: " << result.iterations << std::endl;

    // 验证物质守恒
    std::vector<double> totalMoles(3, 0.0);
    for (const auto& phase : result.phases) {
        for (size_t i = 0; i < 3; ++i) {
            totalMoles[i] += phase.state.moleNumbers[i];
        }
    }

    for (size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(totalMoles[i], input.feedMoles[i], 1e-6)
            << "Component " << i << " mass balance failed";
    }
}

// Test 3: 相数限制测试
TEST_F(StableMethodTest, PhaseNumberLimit) {
    FlashInput input;
    input.temperature = 250.0;
    input.pressure = 2e6;
    input.feedMoles = {0.5, 0.3, 0.2};

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 10;
    opt.max_phases = 2;  // 限制最多2相
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, elementMatrix, opt);

    EXPECT_TRUE(result.success);
    EXPECT_LE(result.numPhases(), 2) << "应遵守max_phases限制";

    std::cout << "Test 3 - Phase Number Limit: "
              << result.numPhases() << " phases (max=2)" << std::endl;
}

// Test 4: 迭代限制测试
TEST_F(StableMethodTest, IterationLimit) {
    FlashInput input;
    input.temperature = 250.0;
    input.pressure = 2e6;
    input.feedMoles = {0.5, 0.3, 0.2};

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 2;  // 限制最多2次外层迭代
    opt.max_phases = 5;
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, elementMatrix, opt);

    EXPECT_TRUE(result.success);

    std::cout << "Test 4 - Iteration Limit: "
              << result.numPhases() << " phases (max_iter=2)" << std::endl;
}

// Test 5: 与快速法对比
TEST_F(StableMethodTest, CompareWithFastMethod) {
    FlashInput input;
    input.temperature = 250.0;
    input.pressure = 2e6;
    input.feedMoles = {0.5, 0.3, 0.2};

    // 快速法
    SolveOptions optFast;
    optFast.phase_determination_strategy = "fast";
    optFast.enable_stability_test = true;

    auto resultFast = flash->solve(input, elementMatrix, optFast);

    // 稳定法
    SolveOptions optStable;
    optStable.phase_determination_strategy = "stable";
    optStable.enable_stability_test = true;
    optStable.max_stability_iterations = 10;
    optStable.max_phases = 3;

    auto resultStable = flash->solve(input, elementMatrix, optStable);

    EXPECT_TRUE(resultFast.success);
    EXPECT_TRUE(resultStable.success);

    std::cout << "Test 5 - Fast vs Stable: "
              << "Fast=" << resultFast.numPhases() << " phases, "
              << "Stable=" << resultStable.numPhases() << " phases" << std::endl;

    // 两种方法应得到相同或相近的相数
    EXPECT_LE(std::abs(static_cast<int>(resultFast.numPhases()) -
                       static_cast<int>(resultStable.numPhases())), 1)
        << "快速法和稳定法结果差异过大";
}

// Test 6: Gibbs能改善阈值测试
TEST_F(StableMethodTest, GibbsImprovementThreshold) {
    FlashInput input;
    input.temperature = 250.0;
    input.pressure = 2e6;
    input.feedMoles = {0.5, 0.3, 0.2};

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 10;
    opt.max_phases = 5;
    opt.gibbs_improvement_threshold = 1e3;  // 较大的阈值，应提前停止
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, elementMatrix, opt);

    EXPECT_TRUE(result.success);

    std::cout << "Test 6 - Gibbs Threshold: "
              << result.numPhases() << " phases (threshold=1e3)" << std::endl;
}

// Test 7: 低温VLE（更容易分相）
TEST_F(StableMethodTest, LowTemperatureVLE) {
    FlashInput input;
    input.temperature = 200.0;  // K - 低温
    input.pressure = 1e6;       // Pa
    input.feedMoles = {0.4, 0.4, 0.2};

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 10;
    opt.max_phases = 3;
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, elementMatrix, opt);

    EXPECT_TRUE(result.success);
    EXPECT_GE(result.numPhases(), 1);

    std::cout << "Test 7 - Low Temperature VLE: "
              << result.numPhases() << " phases at T=" << input.temperature << "K" << std::endl;
}

// Test 8: 含水三相体系（VLLE）- 验证水体系专用初始化
TEST_F(StableMethodTest, WaterSystemVLLE) {
    // 重新创建水体系的热力学后端
    backend = std::make_unique<thermo::ThermoPackBackend>("H2O,C1,nC6", "PR");
    flash = std::make_unique<RandFlash>(*backend, *linSolver);

    // 单位元素矩阵（非反应体系）
    std::vector<std::vector<double>> waterElementMatrix = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    FlashInput input;
    input.temperature = 298.15;  // K - 常温
    input.pressure = 101325.0;   // Pa - 1 atm
    input.feedMoles = {0.45, 0.05, 0.50};  // 45% H2O, 5% C1, 50% nC6

    SolveOptions opt;
    opt.phase_determination_strategy = "stable";
    opt.enable_stability_test = true;
    opt.max_stability_iterations = 10;
    opt.max_phases = 3;
    opt.verbose_stability_loop = false;

    auto result = flash->solve(input, waterElementMatrix, opt);

    EXPECT_TRUE(result.success);
    EXPECT_EQ(result.numPhases(), 3);  // 期望得到 VLLE（气相 + 油相 + 水相）

    // 验证相分率合理性
    double total_beta = 0.0;
    for (size_t j = 0; j < result.numPhases(); ++j) {
        total_beta += result.beta(j);
    }
    EXPECT_NEAR(total_beta, 1.0, 1e-6);

    std::cout << "Test 8 - Water System VLLE: "
              << result.numPhases() << " phases (expected 3 for VLLE)" << std::endl;

    // 打印各相组成以验证物理合理性
    for (size_t j = 0; j < result.numPhases(); ++j) {
        std::cout << "  Phase " << j << " (beta=" << result.beta(j) << "): ";
        for (size_t i = 0; i < result.phases[j].x.size(); ++i) {
            std::cout << result.phases[j].x[i] << " ";
        }
        std::cout << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

