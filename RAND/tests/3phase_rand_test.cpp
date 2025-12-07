// RAND/tests/3phase_rand_test.cpp
#include <iostream>
#include <vector>
#include <memory>
#include <iomanip>

#include "rand_flash.hpp"

// 简单的辅助打印函数
void printPhase(const std::string& name, const std::vector<double>& n, double totalFeed) {
    double beta = 0.0;
    for (double v : n) beta += v;
    
    std::cout << "  " << std::left << std::setw(10) << name 
              << " Beta: " << std::fixed << std::setprecision(5) << (beta/totalFeed) 
              << " | Comp: ";
    if (beta > 1e-10) {
        for (double v : n) std::cout << std::setprecision(4) << (v/beta) << " ";
    } else {
        std::cout << "(Disappeared)";
    }
    std::cout << "\n";
}

int main() {
  try {
    // 1) 构建 ThermoPack 后端
    // 典型的 LLV 体系：水 + 甲烷 + 正己烷
    // H2O: 水相主体
    // C1: 气相主体
    // nC6: 油相主体
    std::string comps = "H2O,C1,nC6";

    // 注意：ThermoPack 通常支持 "H2O" 这种写法，如果报错请尝试 "WATER"
    thermo::ThermoPackBackend backend(comps, "SRK", "vdW", "Classic", "Default", false);

    auto linSolverPtr = ls::createEigenSolver();
    randflash::RandFlash flash(backend, *linSolverPtr);

    // 2) 进料条件：298.15 K, 1 atm
    double T = 298.15;     
    double P = 1.01325e5;  // Pa
    
    // 进料组成：水 0.45, 甲烷 0.05, 正己烷 0.5
    // 这样应该形成：一层水，一层油，和少量的气
    std::vector<double> feed = { 0.45, 0.05, 0.5 };
    
    // 归一化进料（防守）
    double sum_feed = 0.0;
    for (double v : feed) sum_feed += v;
    for (double& v : feed) v /= sum_feed;

    randflash::FlashInput input{T , P , feed};

    // 3) ElementMatrix (单位阵，无反应)
    size_t C = feed.size();
    std::vector<std::vector<double>> elementMatrix(C, std::vector<double>(C, 0.0));
    for (size_t i = 0; i < C; ++i) elementMatrix[i][i] = 1.0;

    // 4) 【关键】手动构造三相初值 (避免相同初值导致矩阵奇异)
    // 组分顺序: 0:H2O, 1:C1, 2:nC6
    
    // 猜测相 1 (Vapor): 富含 C1
    std::vector<double> guessV = { 0.02, 0.78, 0.2 }; 
    
    // 猜测相 2 (Liquid 1 - Oil): 富含 nC6
    std::vector<double> guessL1 = { 0.02, 0.0004, 0.9796 };
    
    // 猜测相 3 (Liquid 2 - Aqueous): 富含 H2O
    std::vector<double> guessL2 = { 0.9999, 0.00005, 0.00005 };

    std::vector<std::vector<double>> initialComps = { guessV, guessL1, guessL2 };

    std::cout << "Starting 3-Phase Flash Calculation (Water + Methane + n-Hexane)...\n";

    // 5) 调用多相闪蒸接口
    int    maxIter = 20;
    double tol     = 1e-8;

    auto res = flash.solveMultiPhase(
      input,
      elementMatrix,
      initialComps, // 传入显式初值
      maxIter,
      tol
    );

    if (!res.success) {
      std::cerr << "[Test] 3-Phase Flash FAILED to converge.\n";
      // 即使失败也打印一下最后结果用于调试
    } else {
        flash.printResult(res);
    }

  } catch (const std::exception& e) {
    std::cerr << "[Test] Exception: " << e.what() << "\n";
    return 1;
  }
}