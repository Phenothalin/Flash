// tests/randflash_multi_systems_test.cpp
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <string>
#include <numeric>
#include <iostream>

#include "rand_flash.hpp"
// #include "thermopack_backend.hpp" // 包含你的实际后端头文件

class MultiSystemFlashTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;

    void SetUp() override {
        // 只初始化通用的线性求解器
        // 具体的 ThermoPack 后端因为组分不同，需要在每个测试中单独构建
        linSolver = ls::createEigenSolver();
    }

    void TearDown() override {
        // 自动清理
    }

    // 辅助函数：归一化
    void normalize(std::vector<double>& v) {
        double sum = std::accumulate(v.begin(), v.end(), 0.0);
        if (sum > 1e-12) {
            for (auto& val : v) val /= sum;
        }
    }

    /**
     * @brief 核心测试函数：接收组分列表和条件，构建新的后端并运行闪蒸
     * * @param compStr  组分字符串，例如 "C1,C2"
     * @param T        温度 (K)
     * @param P        压力 (Pa)
     * @param z        进料组成
     * @param caseName 测试案例名称
     */
    void RunSystemTest(const std::string& compStr, 
                       double T, 
                       double P, 
                       std::vector<double> z, 
                       const std::string& caseName) {
        
        // 1. 针对当前组分构建后端
        // 注意：每次测试重新构建后端会有一定开销，但为了测试不同物质是必须的
        auto backend = std::make_unique<thermo::ThermoPackBackend>(
            compStr, "SRK", "vdW", "Classic", "Default", false);

        // 2. 构建 Flash 计算器
        randflash::RandFlash flash(*backend, *linSolver);

        // 3. 准备输入
        normalize(z);
        randflash::FlashInput input{T, P, z};
        size_t C = z.size();

        // 构造单位矩阵
        std::vector<std::vector<double>> elementMatrix(C, std::vector<double>(C, 0.0));
        for (size_t i = 0; i < C; ++i) elementMatrix[i][i] = 1.0;

        // 4. 执行计算
        int maxIter = 30;
        double tol = 1e-8;

        auto res = flash.solveMultiPhase(input, elementMatrix, {}, maxIter, tol);
        if(res.success){
            flash.printResult(res);
        }
        // 5. 验证
        EXPECT_TRUE(res.success) << "[" << caseName << "] Flash calculation failed.";
    }
};

// ==========================================
// 不同物质体系的测试用例
// ==========================================

// 案例 1: 经典的二元混合物 (甲烷 + 正癸烷)
// 这是一个性质差异很大的体系，常用于测试 VLE
TEST_F(MultiSystemFlashTest, Binary_C1_nC10) {
    std::string comps = "C1,nC10";
    
    // 条件：300K, 50 bar
    // 甲烷是气体，癸烷是液体，这个条件应该处于气液共存区
    double T = 300.0;
    double P = 5.0e6; 
    
    std::vector<double> z = {0.6, 0.4}; // 60% C1, 40% nC10

    RunSystemTest(comps, T, P, z, "Binary_C1_nC10");
}

// 案例 2: 简单 LPG (液化石油气) 体系 (丙烷 + 正丁烷)
// 性质相近的二元体系
TEST_F(MultiSystemFlashTest, LPG_C3_nC4) {
    std::string comps = "C3,nC4";
    
    // 条件：300K, 5 bar (0.5 MPa)
    // 丙烷和丁烷在这个压力下很容易液化，但也容易气化
    double T = 300.0;
    double P = 0.5e6; 
    
    std::vector<double> z = {0.5, 0.5};

    RunSystemTest(comps, T, P, z, "LPG_C3_nC4");
}

// 案例 3: 含非烃类气体体系 (CO2 + N2)
// 低温测试，这是典型的气体分离工况
TEST_F(MultiSystemFlashTest, GasMix_CO2_N2) {
    std::string comps = "N2,CO2";
    
    // 条件：230K (-43C), 20 bar
    // CO2 在这个温度下接近液态/固态边界，N2 是气态
    double T = 230.0;
    double P = 2.0e6; 
    
    std::vector<double> z = {0.2, 0.8}; 

    RunSystemTest(comps, T, P, z, "Cryo_N2_CO2");
}

// 案例 4: 酸气/酸性气体体系 (H2S + C1 + CO2)
// 模拟天然气开采中的高酸气工况
TEST_F(MultiSystemFlashTest, SourGas_H2S_Mix) {
    std::string comps = "C1,CO2,H2S";
    
    // 条件：280K, 40 bar
    double T = 280.0;
    double P = 4.0e6;
    
    std::vector<double> z = {0.5, 0.2, 0.3};

    RunSystemTest(comps, T, P, z, "SourGas_Mix");
}

// 案例 5: 三元烃类体系 (C1 + C3 + nC7)
// 轻、中、重三个组分
TEST_F(MultiSystemFlashTest, Ternary_Hydrocarbon) {
    std::string comps = "C1,C3,nC7";
    
    double T = 320.0;
    double P = 3.0e6; // 30 bar
    
    // 少量庚烷，大量甲烷丙烷
    std::vector<double> z = {0.7, 0.2, 0.1};

    RunSystemTest(comps, T, P, z, "Ternary_C1_C3_C7");
}

// 案例 6: 原始的复杂 11 组分体系 (回归测试)
// 确保新架构也能跑通原来的复杂例子
TEST_F(MultiSystemFlashTest, Original_11_Comps) {
    std::string comps = "N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7";
    
    double T = 295.0;
    double P = 2.0e6;
    
    std::vector<double> z = {
      1, 0.015, 0.55, 0.14, 0.12,
      0.05, 0.045, 0.03, 0.025, 0.012, 0.01
    };

    RunSystemTest(comps, T, P, z, "Original_Complex_Mix");
}

//LLE TEST
TEST_F(MultiSystemFlashTest, LLE_Water_Hexane) {

    // 典型的液液分层体系：水 + 正己烷
    // H2O (极性) 和 nC6 (非极性) 在常温下互不相溶
    std::string comps = "H2O,nC6";

    // 温度 298.15 K (25°C)
    double T = 298.15;

    // 压力 1 atm (1.01325 bar)，确保在液相区
    double P = 1.01325e5;

    // 摩尔组分：50% 水，50% 正己烷
    // 这种配比将导致明显的两相分裂（富水相 和 富油相）
    std::vector<double> z = {
        0.5, 0.5
    };

    RunSystemTest(comps, T, P, z, "LLE_Water_Hexane_Mix");
}

// main 函数
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}