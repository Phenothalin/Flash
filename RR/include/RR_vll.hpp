#pragma once
#include <vector>
#include <stdexcept>
#include "thermo_backend.hpp"

namespace rr_vll {

// 三相VLL闪蒸结果结构体
struct VLLResult {
  // 相分率
  double vapor_fraction;      // β_V
  double liquid1_fraction;    // β_L1
  double liquid2_fraction;    // β_L2 (参考相)

  // 相组成
  std::vector<double> vapor_comp;    // y_i
  std::vector<double> liquid1_comp;  // x_L1,i
  std::vector<double> liquid2_comp;  // x_L2,i (参考相)

  // 收敛信息
  bool converged;
  int iterations;
  double convergence_error;

  VLLResult()
    : vapor_fraction(0.0), liquid1_fraction(0.0), liquid2_fraction(0.0),
      converged(false), iterations(0), convergence_error(0.0) {}
};

// 三相VLL闪蒸求解器
class VLLFlash {
public:
  // 构造函数
  // P, T: 压力和温度
  // z: 总组成
  // thermo_backend: 热力学后端
  // initial_vapor_comp, initial_liquid1_comp, initial_liquid2_comp: 初始相组成猜测
  VLLFlash(double pressure,
           double temperature,
           const std::vector<double>& z,
           thermo::IThermoBackend& thermo_backend,
           const std::vector<double>& initial_vapor_comp = {},
           const std::vector<double>& initial_liquid1_comp = {},
           const std::vector<double>& initial_liquid2_comp = {});

  // 执行三相闪蒸计算
  VLLResult calculate(double tolerance = 1e-6, int max_iterations = 100);

  // Print formatted 3-phase flash results
  void printResult(const VLLResult& result) const;

private:
  // 输入参数
  double P_;
  double T_;
  std::vector<double> z_;
  thermo::IThermoBackend& thermo_backend_;

  // 当前状态
  std::vector<double> y_;      // 气相组成
  std::vector<double> x_L1_;   // 液相1组成
  std::vector<double> x_L2_;   // 液相2组成 (参考相)
  double beta_V_;              // 气相分率
  double beta_L1_;             // 液相1分率

  // K值 (相对于L2参考相)
  std::vector<double> K_V_;    // K_V,i = y_i / x_L2,i
  std::vector<double> K_L1_;   // K_L1,i = x_L1,i / x_L2,i

  // 求解耦合的RR方程组
  bool solveRRSystem(double tolerance, int max_iterations);

  // 从β和K计算相组成
  void computePhaseCompositions();

  // 从逸度系数更新K值
  void updateKValues();

  // 检查收敛性
  bool checkConvergence(const std::vector<double>& K_V_old,
                        const std::vector<double>& K_L1_old,
                        double beta_V_old,
                        double beta_L1_old,
                        double tolerance);

  // 初始化K值和相分率
  void initialize();
};

} // namespace rr_vll
