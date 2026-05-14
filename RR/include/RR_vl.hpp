#pragma once
#include <vector>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <limits>
#include "thermo_backend.hpp"

// 闪蒸类型枚举
enum class FlashType { PT, PV, TV };

// 收敛方法枚举
enum class ConvergenceMethod { NEWTON_RAPHSON, HALLEY };

// 牛顿-拉夫逊求解器
class NewtonRaphsonSolver {
public:
  static double solve(const std::function<double(double)> &func,
                      const std::function<double(double)> &func_deriv,
                      double initial_guess, double tol, int max_iter,
                      std::vector<double> *residual_history = nullptr);
};

// 基类 Flash
class Flash {
protected:
  double pressure_;
  double temperature_;
  double vapor_fraction_;
  double initial_T_;
  double initial_P_;
  double initial_V_;
  std::vector<double> composition_;
  std::vector<double> initial_k_;
  std::vector<double> vap_comp_frac_;
  std::vector<double> liq_comp_frac_;
  thermo::IThermoBackend &thermo_backend_;
  FlashType flash_type_;
  int iterations_ = 0;
  std::vector<double> iter_k_history;  // max_k_diff per outer iteration
  std::vector<double> iter_beta_residual_history;  // flattened inner Newton residuals
  std::vector<int> outer_inner_start_indices;      // start offset of each outer iteration in flattened inner history
  std::vector<int> outer_inner_counts;             // number of inner Newton iterations used by each outer iteration

public:
  Flash(const std::vector<double> &composition,
        thermo::IThermoBackend &thermo_backend,
        double pressure, double temperature, double vapor_fraction,
        const std::vector<double> &initial_k = {},
        double initial_T = std::numeric_limits<double>::quiet_NaN(),
        double initial_P = std::numeric_limits<double>::quiet_NaN());

  void validateInput();
  void initializeKWithWilson(double temperature, double pressure);
  void initializeTWithWilson(double pressure);
  void initializePWithWilson(double temperature);

  std::vector<double> computeVapCompFractions(double v_new,
                                              const std::vector<double> &k_current);
  std::vector<double> computeLiqCompFractions(double v_new,
                                              const std::vector<double> &k_current);
  std::vector<double> updateKValues(const std::vector<double> &phi_liquid,
                                    const std::vector<double> &phi_vapor);
  void displayResults() const;

  // Print formatted 2-phase flash results
  void printResult() const;

  // 多相扩展接口（预留）
  virtual void proposeNewPhaseByStabilityAnalysis();

  virtual void calculate(ConvergenceMethod method) = 0;

  // Getters
  [[nodiscard]] auto getVaporFraction() const -> double {
    return vapor_fraction_;
  }
  [[nodiscard]] auto getTemperature() const -> double { return temperature_; }
  [[nodiscard]] auto getPressure() const -> double { return pressure_; }
  [[nodiscard]] auto getIterations() const -> int { return iterations_; }
  [[nodiscard]] auto getIterKHistory() const -> const std::vector<double>& { return iter_k_history; }
  [[nodiscard]] auto getIterBetaResidualHistory() const -> const std::vector<double>& { return iter_beta_residual_history; }
  [[nodiscard]] auto getOuterInnerStartIndices() const -> const std::vector<int>& { return outer_inner_start_indices; }
  [[nodiscard]] auto getOuterInnerCounts() const -> const std::vector<int>& { return outer_inner_counts; }
  auto getVapComp() -> std::vector<double> { return vap_comp_frac_; }
  auto getLiqComp() -> std::vector<double> { return liq_comp_frac_; }
};

// PTFlash 类
class PTFlash : public Flash {
public:
  PTFlash(double pressure, double temperature, const std::vector<double> &composition,
          thermo::IThermoBackend &thermo_backend,
          const std::vector<double> &initial_k = {});
  void calculate(ConvergenceMethod method) override;

private:
  void calculateVaporFraction(ConvergenceMethod method);
  double calculateRachfordRice(const double &vapfrac, const std::vector<double> &K_CURRENT);
  double calculateRachfordRiceDeriv(const double &vapfrac, const std::vector<double> &K_CURRENT);
  double calculateRachfordRiceSecondDeriv(const double &vapfrac, const std::vector<double> &K_CURRENT);
  bool checkConvergence(const std::vector<double> &k_current,
                        const std::vector<double> &k_previous,
                        double v_new, double v_current,
                        double TOL_OUTER, int outer_iter);
};

// PVFlash 类
class PVFlash : public Flash {
public:
  PVFlash(double pressure, double vapor_fraction, const std::vector<double> &composition,
          thermo::IThermoBackend &thermo_backend,
          const double &initial_T = std::numeric_limits<double>::quiet_NaN(),
          const std::vector<double> &initial_k = {});
  void calculate(ConvergenceMethod method) override;

private:
  void calculateTemperature(ConvergenceMethod method);
  double calculateRachfordRice(const std::vector<double> &K_value,
                               const double &vaporfraction,
                               const std::vector<double> &composition);
  double calculateRachfordRiceDeriv(const double &temperature,
                                    const std::vector<double> &K_value);
  double updateTemperature(const double &temperature, const double &err, const double &derr);
  bool checkConvergence(const std::vector<double> &k_current,
                        const std::vector<double> &k_previous,
                        const double &T_current, const double &T_previous,
                        int outer_iter);
};

// TVFlash 类
class TVFlash : public Flash {
public:
  TVFlash(double temperature, double vapor_fraction, const std::vector<double> &composition,
          thermo::IThermoBackend &thermo_backend,
          const double &initial_P = std::numeric_limits<double>::quiet_NaN(),
          const std::vector<double> &initial_k = {});
  void calculate(ConvergenceMethod method) override;

private:
  void calculatePressure(ConvergenceMethod method);
  double calculateRachfordRice(const std::vector<double> &K_value,
                               const double &vaporfraction,
                               const std::vector<double> &composition);
  double calculateRachfordRiceDeriv(const double &pressure,
                                    const std::vector<double> &K_value);
  double updatePressure(const double &pressure, const double &err, const double &derr);
  bool checkConvergence(const std::vector<double> &k_current,
                        const std::vector<double> &k_previous,
                        const double &P_current, const double &P_previous,
                        int outer_iter);
};
