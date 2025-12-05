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
  std::vector<std::vector<double>> n_phase;            // F x C
  std::vector<double> beta;                            // F
  int iterations = 0;                                  
  double mu_infinity_norm = 0.0;                       
  double elem_residual_inf = 0.0;                      
};

class RandFlash {
public:
  RandFlash(thermo::IThermoBackend& thermo,
    ls::LinearSolverInterface& linearSolver);

  // 保持旧接口用于兼容性测试
  FlashResult solveTwoPhase(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& initialVaporComposition = {},
    const std::vector<double>& initialLiquidComposition = {},
    int maxIterations = 50,
    double tolerance = 1e-8);
      
  MultiFlashResult solveMultiPhaseCore(
    const FlashInput& input,
    int nPhases,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions = {},
    int maxIterations = 50,
    double tolerance = 1e-8);

  MultiFlashResult solveMultiPhase(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions = {},
    int maxIterations = 50,
    double tolerance = 1e-8);  
      
  struct ConvergenceInfo {
    double max_mu_diff;   
    double elem_error;    
    bool   converged;     
  };

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
  double applyUpdate(
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
  // 所有的 Newton 迭代逻辑移到这里，solveTwoPhase 和 solveMultiPhaseCore 只是它的包装
  MultiFlashResult solveGeneral(
    SystemContext& sys,
    int maxIterations,
    double tolerance);
};

} // namespace randflash