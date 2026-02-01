// RAND/include/rand_flash.hpp

#pragma once
#include "thermo_backend.hpp"
#include "linear_solver.hpp"
#include <vector>
#include <memory>

namespace randflash {

// 通用初始化结果（不再限于两相）
struct InitResult {
  std::vector<std::vector<double>> n_phases; // [F][C]
  std::vector<double> beta;                  // [F]
  std::vector<std::vector<double>> compositions; // [F][C] (x, y...)
};

struct FlashInput {
  double temperature = 0.0;
  double pressure    = 0.0;
  std::vector<double> feedMoles;  // size C
};

// 控制 solve() 行为的选项。
//
// 设计目标：
// 1) 默认保持现有行为：enable_stability_test=true 且不提供初值时自动做 stability analysis。
// 2) 允许关闭 stability analysis，并强制指定相数（reaction 体系常见：固定相数后在元素守恒约束下做 Gibbs 最小化）。
// 3) 允许用户给定初始相组成（此时无论 enable_stability_test 与否，均跳过 stability analysis）。
struct SolveOptions {
  // 是否启用相稳定性分析。
  // - true ：若未提供 initial_phase_compositions，则走原自动流程。
  // - false：若未提供 initial_phase_compositions，则必须提供 forced_phase_count。
  bool enable_stability_test = true;

  // 当 enable_stability_test=false 且 initial_phase_compositions 为空时，强制相数。
  // 取值：F>=1。
  int forced_phase_count = 0;

  // 强制相型（phase flag）。若为空，将默认采用：
  // - F==1：liquid
  // - F>=2：第 0 相 vapor，其余 liquid
  std::vector<int> forced_phase_flags;

  // 初始相组成猜测（每行一个相的组成向量，需归一化或近似归一化）。
  // 非空时：直接按给定相数求解并跳过 stability analysis。
  std::vector<std::vector<double>> initial_phase_compositions;
};

struct PhaseFixOptions {
  double eig_floor     = 1e-10;
  double eps_shift     = 1e-12;
  double inv_tol       = 1e-12;
  int    max_repair_it = 3;
};

struct PhaseFixResult {
  bool   applied = false;
  double lam_min_before = 0.0;
  double lam_min_after  = 0.0;
  std::vector<std::vector<double>> m_fixed;
  std::vector<std::vector<double>> M_fixed;
};

struct PhaseContext {
  thermo::PhaseState state;
  std::vector<double> x;                    
  std::vector<std::vector<double>> m;       
  std::vector<std::vector<double>> M;       
  std::vector<double> mu;                   
};

struct SystemContext {
  double temperature = 0.0;
  double pressure    = 0.0;
  std::vector<double> feedMoles;                    
  std::vector<std::vector<double>> elementMatrix;   
  std::vector<PhaseContext> phases;    // size F
};

// 保持 FlashResult 兼容旧的 solveTwoPhase 接口
struct FlashResult {
  bool success = false;
  double pressure = 0.0;
  double temperature = 0.0;
  std::vector<double> feedComposition;
  std::vector<double> vaporComposition;
  std::vector<double> liquidComposition;
  double vaporFraction = 0.0;
  size_t iterations = 0;
  double convergenceError = 0.0;
};

struct MultiFlashResult {
  bool success = false;
  double pressure = 0.0;
  double temperature = 0.0;
  std::vector<double> feedComposition;
  std::vector<PhaseContext> phases;                    // F - 包含各相的完整上下文（含phaseFlag、moleNumbers等）
  int iterations = 0;
  double mu_infinity_norm = 0.0;
  double elem_residual_inf = 0.0;

  // 便捷方法：获取相数
  size_t numPhases() const { return phases.size(); }

  // 便捷方法：获取第j相的总摩尔数（beta）
  double beta(size_t j) const {
    if (j >= phases.size()) return 0.0;
    double sum = 0.0;
    for (double n : phases[j].state.moleNumbers) sum += n;
    return sum;
  }

  // 便捷方法：获取第j相的摩尔数向量
  const std::vector<double>& n_phase(size_t j) const {
    static const std::vector<double> empty;
    if (j >= phases.size()) return empty;
    return phases[j].state.moleNumbers;
  }
};

class RandFlash {
public:
  RandFlash(thermo::IThermoBackend& thermo,
    ls::LinearSolverInterface& linearSolver);

  // === 通用接口：自动相稳定性判别 + 初始化 + 相分裂 ===
  // 若 initialPhaseCompositions 非空，则跳过相稳定性分析，按用户给定相数求解。
  MultiFlashResult solve(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions = {},
    int maxIterations = 50,
    double tolerance = 1e-8);

  // 新接口：通过 SolveOptions 控制是否执行 stability analysis / 是否强制相数。
  // 说明：
  // - opt.initial_phase_compositions 非空时：直接按给定相数求解（跳过 stability）。
  // - opt.enable_stability_test=true 且 initial 为空：保持原自动流程。
  // - opt.enable_stability_test=false 且 initial 为空：必须给 forced_phase_count。
  MultiFlashResult solve(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations = 50,
    double tolerance = 1e-8);

  // 保持旧接口用于兼容性测试
  FlashResult solveTwoPhase(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& initialVaporComposition = {},
    const std::vector<double>& initialLiquidComposition = {},
    int maxIterations = 50,
    double tolerance = 1e-8);

  // 注意：v5 起不再暴露 solveMultiPhase / solveMultiPhaseCore。
  // 统一通过 solve() 使用：
  // - 不提供 initialPhaseCompositions：内部做相稳定性分析 + 自动选相数/相型
  // - 提供 initialPhaseCompositions：按给定相数组合做初始化并求解
      
  struct ConvergenceInfo {
    double max_mu_diff;   
    double elem_error;    
    bool   converged;     
  };

  void printResult(const MultiFlashResult& res) const;

private:
  thermo::IThermoBackend& thermo_;
  ls::LinearSolverInterface& linearSolver_;

  void updatePhaseChemistry(PhaseContext& phaseCtx);

  void fixPhaseHessian(PhaseContext& phaseCtx,
                       const PhaseFixOptions& opt);

  void assembleLocalJacobian(const thermo::PhaseState& state,
    std::vector<std::vector<double>>& m,
    std::vector<double>& mu);

  // 全局系统组装（已支持多相，无需大改，只需确认内部逻辑）
  void assembleGlobalSystem(
    double temperature,
    const std::vector<std::vector<std::vector<double>>>& Ms,
    const std::vector<std::vector<double>>& mus,
    const std::vector<std::vector<double>>& nPhases,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& feedComposition, 
    std::vector<double>& Acoef,
    std::vector<double>& rhs); 
    
  // 两相初始化保持原样，仅供 solveTwoPhase 调用
  InitResult initializeTwoPhase(
    const SystemContext& sys,                
    const std::vector<double>& vaporGuess,                 
    const std::vector<double>& liquidGuess,                
    double margin = 1e-3) const;     

  InitResult initializeThreePhaseWater(
    const SystemContext& sys, 
    double margin = 1e-10) const;

  // 通用三相初始化：不假设必须含水
  InitResult initializeThreePhaseGeneric(
    const SystemContext& sys,
    double margin = 1e-10) const;

  // 用"相组成猜测"构造满足守恒的初始 n（最小二乘 beta + 行缩放严格守恒）
  InitResult initializeFromCompositions(
    const SystemContext& sys,
    int nPhases,
    const std::vector<std::vector<double>>& phaseCompositions,
    const std::vector<int>& phaseFlags,
    double min_phase_moles_ratio = 1e-8) const;

  // ========== Reactive System Initialization ==========

  // 单相反应体系初始化：从元素摩尔数构造可行组成
  // 使用最小二乘法求解 A·n = b，确保元素守恒
  InitResult initializeReactiveSinglePhase(
    const SystemContext& sys,
    const std::vector<double>& elementMoles,
    int phaseFlag) const;

  // 多相反应体系初始化：分配元素到各相
  // 策略：按相型分配初始组成，确保总元素守恒
  InitResult initializeReactiveMultiPhase(
    const SystemContext& sys,
    int numPhases,
    const std::vector<int>& phaseFlags,
    const std::vector<double>& elementMoles) const;

  // === 【修改 1】通用化 LineSearch ===
  double lineSearch(
    const std::vector<std::vector<double>>& nPhases,
    const std::vector<std::vector<double>>& dnPhases,
    const std::vector<std::vector<double>>& gPhases,   
    double init_alpha,             
    double min_alpha,              
    double shrink
  );

  std::vector<double> solveGlobalLinearSystem(
    const std::vector<double>& Acoef,
    const std::vector<double>& rhs,
    double* residual_out = nullptr) const; 
  
  void backSubstituteDeltas(
    double temperature,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<std::vector<double>>>& Ms,
    const std::vector<std::vector<double>>& mus,
    const std::vector<std::vector<double>>& nPhases,
    const std::vector<double>& Lambda,         
    const std::vector<double>& deltaBeta,      
    std::vector<std::vector<double>>& dnPhases 
  ) const;

  // === 【修改 2】通用化 applyUpdate ===
  // 返回 alpha
  void applyUpdate(
    double temperature,
    const std::vector<std::vector<double>>& mus,
    const std::vector<std::vector<double>>& dnPhases,
    std::vector<std::vector<double>>& nPhases_inout);
  
  // === 【修改 3】通用化 checkConvergence ===
  ConvergenceInfo checkConvergence(
    int iternumber,
    double tol,
    double temperature,
    const std::vector<std::vector<double>>& mus,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& nPhases,
    const std::vector<double>& feedComposition) const;

  // === 【新增】通用求解内核 ===
  // 所有的 Newton 迭代逻辑移到这里；solve() / solveTwoPhase() 只是包装。
  MultiFlashResult solveGeneral(
    SystemContext& sys,
    int maxIterations,
    double tolerance);

  // === 内部驱动：给定相数组合猜测 + 相型，完成初始化 + solveGeneral + 结果相序稳定化 ===
  MultiFlashResult solveWithPhaseGuesses(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& phaseCompositions,
    const std::vector<int>& phaseFlags,
    int maxIterations,
    double tolerance);
};

} // namespace randflash