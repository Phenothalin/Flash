#pragma once
#include "thermo_backend.hpp"
#include "linear_solver.hpp"
#include <vector>
#include <memory>

namespace randflash {

struct InitResult {
  std::vector<double> nV, nL;  // 各组分相内摩尔数
  double betaV = 0.0, betaL = 0.0;
  std::vector<double> y0, x0;  // 归一化后的初始 y/x（便于日志）
};

struct PhaseFixOptions {
  double eig_floor     = 1e-10;   // 切空间最小特征值下限 δ
  double eps_shift     = 1e-12;   // 额外小移位 ε
  double inv_tol       = 1e-12;   // 求逆时的数值容差
  int    max_repair_it = 3;       // 失败时再微调次数
};

struct PhaseFixResult {
  bool   applied = false;
  double lam_min_before = 0.0;
  double lam_min_after  = 0.0;
  std::vector<std::vector<double>> m_fixed;
  std::vector<std::vector<double>> M_fixed;
};

struct PhaseContext {
  thermo::PhaseState state;              // T, P, n, phaseFlag
  std::vector<std::vector<double>> m;    // 局部 m 矩阵
  std::vector<double> mu;               // 化学势
};

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
  std::vector<double> feedComposition;                 // size C, absolute moles
  std::vector<std::vector<double>> n_phase;            // F x C, moles in each phase
  std::vector<double> beta;                            // size F, phase totals
  int iterations = 0;                                  // inner Newton iterations (last)
  double mu_infinity_norm = 0.0;                       // max phase-to-phase μ spread
  double elem_residual_inf = 0.0;                      // ||A*(sum n - feed)||_inf
};

class RandFlash {
  public:
    RandFlash(thermo::IThermoBackend& thermo,
      ls::LinearSolverInterface& linearSolver);

    FlashResult solveTwoPhase(
        const thermo::PhaseState& state,
        const std::vector<std::vector<double>>& elementMatrix,
        const std::vector<double>& initialVaporComposition = {},
        const std::vector<double>& initialLiquidComposition = {},
        int maxIterations = 50,
        double tolerance = 1e-8);
        
      // 用于“收敛判断”的返回
      struct ConvergenceInfo {
        double max_mu_diff;   // max |muV[i]-muL[i]|
        double elem_error;    // || A*(nV+nL-feed) ||_inf
        bool   converged;     // (max_mu_diff < tol && elem_error < 1e-8)
      };

      // RandFlash(const RandFlash&) = delete;
      // RandFlash& operator=(const RandFlash&) = delete;
    
      // // 如有需要，允许移动（方便以后放在容器里）
      // RandFlash(RandFlash&&) = default;
      // RandFlash& operator=(RandFlash&&) = default;

  private:
  // std::shared_ptr<thermo::IThermoBackend> thermo_;
  thermo::IThermoBackend& thermo_;
  ls::LinearSolverInterface& linearSolver_;

  void assembleLocalJacobian(const thermo::PhaseState& state,
    std::vector<std::vector<double>>& m,
    std::vector<double>& mu);

  void assembleGlobalSystem(
    double temperature,
    const std::vector<std::vector<double>>& MV,
    const std::vector<std::vector<double>>& ML,
    const std::vector<double>& muV,
    const std::vector<double>& muL,
    const std::vector<double>& moleNumbersV,
    const std::vector<double>& moleNumbersL,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& feedComposition, 
    std::vector<double>& Acoef,
    std::vector<double>& rhs);
  
  InitResult initializeTwoPhase(
    const std::vector<double>& feed,                       // 进料各组分摩尔数
    const std::vector<double>& vaporGuess,                 // 允许空
    const std::vector<double>& liquidGuess,                // 允许空
    double margin = 1e-3) const;                           // 远离边界的小裕度

  double armijoLineSearch(
    const std::vector<double>& nV,
    const std::vector<double>& nL,
    const std::vector<double>& dnV,
    const std::vector<double>& dnL,
    double temperature, double pressure,
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& feed,
    double eta,           // e.g. 1e-2
    double alpha0,        // 1.0
    double c,             // 1e-4
    double shrink
  );
  double lineSearch(
    const std::vector<double>& nV,
    const std::vector<double>& nL,
    const std::vector<double>& dnV,
    const std::vector<double>& dnL,
    const std::vector<double>& gV,   
    const std::vector<double>& gL,   
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
    const std::vector<std::vector<double>>& MV,
    const std::vector<std::vector<double>>& ML,
    const std::vector<double>& muV,
    const std::vector<double>& muL,
    const std::vector<double>& nV,
    const std::vector<double>& nL,
    const std::vector<double>& Lambda,         // 长度 E
    const std::vector<double>& deltaBeta,      // 长度 2: {Δβ_v, Δβ_l}
    std::vector<double>& dnV,
    std::vector<double>& dnL) const;

  double applyUpdate(
    double temperature,
    const std::vector<double>& muV,
    const std::vector<double>& muL,
    const std::vector<double>& dnV,
    const std::vector<double>& dnL,
    std::vector<double>& nV_inout,
    std::vector<double>& nL_inout) ;
    
  ConvergenceInfo checkConvergence(
    int iternumber,
    double tol,
    const std::vector<double>& muV,
    const std::vector<double>& muL,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& nV,
    const std::vector<double>& nL,
    const std::vector<double>& feedComposition) const;
  };

} // namespace randflash
