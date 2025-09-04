#pragma once
#include "thermo_adapter.hpp"
#include "linear_solver.hpp"
#include <vector>
#include <memory>

namespace randflash {

struct InitResult {
  std::vector<double> nV, nL;  // 各组分相内摩尔数
  double betaV = 0.0, betaL = 0.0;
  std::vector<double> y0, x0;  // 归一化后的初始 y/x（便于日志）
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
  RandFlash(
    PropertyPackageType packageType,
    std::shared_ptr<material_object::Cluster> componentCluster,
    ls::LinearSolverInterface& linearSolver);

  FlashResult solveTwoPhase(
    double pressure,
    double temperature,
    const std::vector<double>& feedComposition,
    const std::vector<std::vector<double>>& elementMatrix, // E×C 元素映射
    const std::vector<double>& initialVaporComposition = {},  // 新增：气相初始组成
    const std::vector<double>& initialLiquidComposition = {},// 新增：液相初始组成
    int maxIterations = 50,
    double tolerance    = 1e-8);

private:
  thermo::PropertyPackageAdapter vaporModel_, liquidModel_;
  ls::LinearSolverInterface& linearSolver_;

  void assembleLocalJacobian(
    double temperature,
    double pressure,
    std::vector<double>& moleNumbers,
    thermo::PropertyPackageAdapter& model,
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
    double pressure,
    double temperature,
    const std::vector<double>& feed,                       // 进料各组分摩尔数
    const std::vector<std::vector<double>>& elementMatrix, // A 矩阵
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
};

class RandMultiphase {
  public:
  struct Options {
    int max_newton_iter; // Newton iters for a fixed phase count
    int max_addremove_cycles; // Outer cycles for add/remove phases
    double tol_mu; // μ spread tolerance
    double tol_elem; // element residual tolerance
    double remove_phase_tol_rel; // remove phase if beta < rel * Ftot
    double line_alpha_min; // min line-search step
    double line_shrink; // backtracking shrink
    double jitter_scale; // small diagonal regularization
    int max_phases; // safety cap on phases
    bool verbose; // print progress to std::cout
    
    
    Options()
    : max_newton_iter(50)
    , max_addremove_cycles(6)
    , tol_mu(1e-8)
    , tol_elem(1e-8)
    , remove_phase_tol_rel(1e-12)
    , line_alpha_min(1e-12)
    , line_shrink(0.5)
    , jitter_scale(1e-6)
    , max_phases(5)
    , verbose(true) {}
    };
  
  
  RandMultiphase(
  PropertyPackageType packageType,
  std::shared_ptr<material_object::Cluster> componentCluster,
  ls::LinearSolverInterface& linearSolver,
  const Options& opt = Options{});
  
  
  ~RandMultiphase();
  
  
  // Solve isothermal-isobaric flash with automatic phase count discovery.
  // feed: absolute moles (size C). elementMatrix: E x C.
  // initial_x: optional initial phase compositions (each size C, sum=1). If empty,
  // solver starts from a symmetric 2-phase split with x=z.
  MultiFlashResult solvePT(
  double pressure,
  double temperature,
  const std::vector<double>& feed,
  const std::vector<std::vector<double>>& elementMatrix,
  const std::vector<std::vector<double>>& initial_x = {});
  
  
  private:
  struct Impl; // PImpl to keep header light & ABI-stable
  std::unique_ptr<Impl> impl_; // defined in rand_multiphase.cpp
  };

} // namespace randflash
