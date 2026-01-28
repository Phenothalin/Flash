// tests/randflash_tp_two_phase.cpp
#include <iostream>
#include <vector>
#include <memory>

#include "rand_flash.hpp"

int main() {
  try {
    // 1) 构建 ThermoPack 后端：这里用 Peng-Robinson + vdW + Classic alpha
    // 组分顺序要和你的 feed z 保持一致
    std::string comps =
      "N2,CO2,C1,C2,C3,"
      "iC4,nC4,iC5,nC5,nC6,nC7";

    thermo::ThermoPackBackend backend(
      comps, "SRK", "vdW", "Classic", "Default", false);

    // 2) 线性求解器（和你以前单相/两相测试一致）
    auto linSolverPtr = ls::createEigenSolver();
    randflash::RandFlash flash(backend, *linSolverPtr);

    // 3) 进料条件：T, P, 总摩尔组成（用你之前那个 11 组分例子）
    double T = 295.0;      // K
    double P = 2.0e6;      // Pa
    std::vector<double> feed = {
      1, 0.015, 0.55, 0.14, 0.12,
      0.05, 0.045, 0.03, 0.025, 0.012, 0.01
    };
    double sum_feed = 0.0;
    for (double v : feed) sum_feed += v;
    const size_t C = feed.size();

    // 设总进料为 1 mol（你后面算法只看比例，这里 1 mol 就够）
    // std::vector<double> nFeed(C);
    // double totalFeed = 1.0;
    // for (size_t i = 0; i < C; ++i) nFeed[i] = z[i] * totalFeed;

    // 4) 构造 FlashInput
    randflash::FlashInput input{T , P , feed};

    // 5) elementMatrix：这里用单位矩阵，表示“每个成分各自守恒”
    std::vector<std::vector<double>> elementMatrix(C, std::vector<double>(C, 0.0));
    for (size_t i = 0; i < C; ++i) elementMatrix[i][i] = 1.0;


    int    maxIter = 20;
    double tol     = 1e-8;
    randflash::SolveOptions options;
    auto res = flash.solve(
      input,
      elementMatrix,
      options,
      maxIter,
      tol
    );

    if (!res.success) {
      std::cerr << "[Test] solveTwoPhase did not converge.\n";
      return 1;
    }

    flash.printResult(res);

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "[Test] Exception: " << e.what() << "\n";
    return 1;
  }
}
