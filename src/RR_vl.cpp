#include "RR_vl.hpp"

#include <algorithm>
#include <cmath>

#include "eos_model.hpp"
#include "property_package.hpp"

// 内部优化的牛顿-拉夫逊求解器
// auto NewtonRaphsonSolver::solve(const std::function<double(double)> &func,
//                                 const std::function<double(double)> &func_deriv,
//                                 double initial_guess, double tol,
//                                 int max_iter) -> double {
//   double x_value = initial_guess;
//   double alpha = 1;
//   for (int iter = 0; iter < max_iter; ++iter)
//   {
//     double fx_value = func(x_value);
//     double dfx_value = func_deriv(x_value);

//     if (dfx_value == 0.0)
//     {
//       throw std::runtime_error("牛顿-拉夫逊方法中遇到零导数。");
//     }

//     double x_new = x_value - alpha * (fx_value / dfx_value);

//     std::cout << "迭代 " << iter + 1 << ": x = " << x_value
//               << ", x_new = " << x_new << ", f(x_new) = " << func(x_new)
//               << "\n";

//     if (std::abs(x_new - x_value) < tol)
//     {
//       std::cout << iter + 1 << " 次迭代后收敛。\n";
//       return x_new;
//     }

//     x_value = x_new;
//   }

//   throw std::runtime_error("牛顿-拉夫逊方法在最大迭代次数内未收敛。");
// }

// 普通牛顿-拉夫逊求解器
auto NewtonRaphsonSolver::solve(const std::function<double(double)> &func,
                                const std::function<double(double)> &func_deriv,
                                double initial_guess, double tol,
                                int max_iter) -> double {
  double x_value = initial_guess;
  double alpha = 1;
    double fx_value = func(x_value);
    double dfx_value = func_deriv(x_value);

    if (dfx_value == 0.0)
    {
      throw std::runtime_error("牛顿-拉夫逊方法中遇到零导数。");
    }

    double x_new = x_value - alpha * (fx_value / dfx_value);
    return x_new;
}

// Flash 基类的构造函数实现 - 增强版
Flash::Flash(const std::vector<double> &composition,
             property_package::PropertyPackage &property_package,
             double pressure, double temperature, double vapor_fraction,
             const std::vector<double> &initial_k, double initial_T,
             double initial_P)
    : pressure_(pressure), temperature_(temperature), initial_T_(initial_T),
      initial_P_(initial_P), composition_(composition), initial_k_(initial_k),
      vapor_fraction_(vapor_fraction), initial_V_(0.5), // 默认初始气相分率为0.5
      vap_comp_frac_(composition.size(), 0),
      liq_comp_frac_(composition.size(), 0),
      property_package_(property_package) {
  // 基本的通用输入验证，现在只验证组成
  validateInput();
}

// 基类通用输入验证
void Flash::validateInput() {
  if (!std::isnan(pressure_) && pressure_ <= 0.0)
  {
    throw std::invalid_argument("压力必须大于零。");
  }
  if (!std::isnan(temperature_) && temperature_ <= 0.0)
  {
    throw std::invalid_argument("温度必须大于零。");
  }
  if (!std::isnan(vapor_fraction_) &&
      (vapor_fraction_ < 0.0 || vapor_fraction_ > 1.0))
  {
    throw std::invalid_argument("气相分率必须在0到1之间");
  }
  if (!std::isnan(initial_T_) && initial_T_ <= 0)
  {
    throw std::invalid_argument("温度必须大于零。");
  }
  if (composition_.empty())
  {
    throw std::invalid_argument("进料摩尔分率不能为空。");
  }
  double sum_composition = 0.0;
  for (auto val : composition_)
  {
    if (val < 0.0)
    {
      throw std::invalid_argument("进料摩尔分率不能包含负值。");
    }
    sum_composition += val;
  }
  if (std::abs(sum_composition - 1.0) > 1e-6)
  {
    throw std::invalid_argument("进料摩尔分率之和必须为1。");
  }
  // 如果 initial_k_ 非空，则做相应检查
  if (!initial_k_.empty())
  {
    if (initial_k_.size() != composition_.size())
    {
      throw std::invalid_argument(
          "初始K值向量必须与进料摩尔分率向量具有相同的大小。");
    }
    for (auto k_val : initial_k_)
    {
      if (k_val <= 0.0)
      {
        throw std::invalid_argument("初始K值必须大于零。");
      }
    }
  }
}

// Wilson 方法计算 K（基类通用方法）
void Flash::initializeKWithWilson(double temperature, double pressure) {
  if (pressure <= 0.0)
  {
    throw std::invalid_argument("压力必须大于零。");
  }

  initial_k_.clear(); // 清空旧的K值
  property_package::EosComponentParameters
      eos_correction_coefficient_parameters =
          property_package_.getCriticalParameters();
  for (size_t index = 0; index < composition_.size(); ++index)
  {
    double critical_temp =
        eos_correction_coefficient_parameters.critical_temperatures[index];
    double critical_pressure =
        eos_correction_coefficient_parameters.critical_pressures[index];
    double acentric_factors =
        eos_correction_coefficient_parameters.acentric_factors[index];
    if (critical_temp <= 0.0 || critical_pressure <= 0.0)
    {
      throw std::invalid_argument("Wilson方法需要临界温度和压力大于零。");
    }
    // Wilson 方法
    double ln_ki =
        std::log(critical_pressure / pressure) +
        5.373 * (1.0 + acentric_factors) * (1.0 - critical_temp / temperature);
    double ki_value = std::exp(ln_ki);
    initial_k_.push_back(ki_value);
  }
}

// 初始化初始温度估计（从PVFlash移至基类）
void Flash::initializeTWithWilson(double pressure) {
  if (pressure <= 0.0)
  {
    throw std::invalid_argument("压力必须大于零。");
  }

  property_package::EosComponentParameters
      eos_correction_coefficient_parameters =
          property_package_.getCriticalParameters();
  double T_LOW = 1E100;
  double T_HIGH = 0.0;
  double T_MAX = 50000.0;
  std::vector<double> T_process(composition_.size(), 0.0);
  for (size_t index = 0; index < composition_.size(); ++index)
  {
    double critical_temp =
        eos_correction_coefficient_parameters.critical_temperatures[index];
    double critical_pressure =
        eos_correction_coefficient_parameters.critical_pressures[index];
    double acentric_factors =
        eos_correction_coefficient_parameters.acentric_factors[index];
    if (critical_temp <= 0.0 || critical_pressure <= 0.0)
    {
      throw std::invalid_argument("Wilson方法需要临界温度和压力大于零。");
    }
    T_process[index] = (5.373 * critical_temp * (acentric_factors + 1)) /
                       (5.373 * (acentric_factors + 1) - log(pressure) +
                        log(critical_pressure));
    // std::cout << T_process[index] << std::endl;
    if (T_process[index] < T_LOW)
    {
      T_LOW = T_process[index];
    }
    if (T_process[index] > T_HIGH)
    {
      T_HIGH = T_process[index];
    }
  }
  // std::cout << T_LOW << "  " << T_HIGH << std::endl;
  if (T_LOW <= 0)
  {
    T_LOW = 1E-12;
  }
  if (T_HIGH <= 0)
  {
    throw std::invalid_argument("温度必须大于零。");
  }
  if (T_HIGH < 0.1 * T_MAX)
  {
    initial_T_ = 0.4 * (T_HIGH + T_LOW);
  }
  else
  {
    initial_T_ = 0;
    for (size_t i = 0; i < composition_.size(); ++i)
    {
      initial_T_ +=
          composition_[i] *
          eos_correction_coefficient_parameters.critical_temperatures[i];
    }
    initial_T_ *= 0.666666;
    if (initial_T_ < T_LOW)
    {
      initial_T_ = T_LOW + 1.0;
    }
  }
}

// 新增：初始化压力的方法
void Flash::initializePWithWilson(double temperature) {
  if (temperature <= 0.0)
  {
    throw std::invalid_argument("温度必须大于零。");
  }

  property_package::EosComponentParameters
      eos_correction_coefficient_parameters =
          property_package_.getCriticalParameters();

  // 使用组分临界压力的加权平均值作为初始估计
  initial_P_ = 0.0;
  double p_bubble = 0;
  double p_dew = 0;
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    double critical_pressure =
        eos_correction_coefficient_parameters.critical_pressures[i];
    double critical_temp =
        eos_correction_coefficient_parameters.critical_temperatures[i];
    double acentric_factors =
        eos_correction_coefficient_parameters.acentric_factors[i];

    if (critical_temp <= 0.0 || critical_pressure <= 0.0)
    {
      throw std::invalid_argument("Wilson方法需要临界温度和压力大于零。");
    }

    p_bubble +=
        composition_[i] * critical_pressure *
        exp(5.373 * (1 + acentric_factors) * (1 - critical_temp / temperature));
    p_dew += composition_[i] /
             (critical_pressure * exp(5.373 * (1 + acentric_factors) *
                                      (1 - critical_temp / temperature)));
    // 基于Wilson公式的逆推，用于估计压力
    //    double term =
    //        5.373 * (1.0 + acentric_factors) * (1.0 - critical_temp /
    //        temperature);
    //    double p_estimate = critical_pressure * std::exp(-term);
    //
    //    initial_P_ += composition_[i] * p_estimate;
  }
  // std::cout<<"p_bubble "<<p_bubble<<"    p_dew "<<p_dew<<"\n";
  initial_P_ = p_bubble + vapor_fraction_ * (p_dew - p_bubble);
  // std::cout<<"initial_p_ "<<initial_P_<<"\n";
  //  确保初始压力为正值
  if (initial_P_ <= 0.0)
  {
    throw std::runtime_error("计算出的初始压力<=0");
  }
}

// 计算气相组成（基类通用方法）
auto Flash::computeVapCompFractions(
    double v_new, const std::vector<double> &k_current) -> std::vector<double> {
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    vap_comp_frac_[i] =
        composition_[i] * k_current[i] / (1.0 + v_new * (k_current[i] - 1.0));
  }
  return vap_comp_frac_;
}

// 计算液相组成（基类通用方法）
auto Flash::computeLiqCompFractions(
    double v_new, const std::vector<double> &k_current) -> std::vector<double> {
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    liq_comp_frac_[i] = composition_[i] / (1.0 + v_new * (k_current[i] - 1.0));
  }
  return liq_comp_frac_;
}

// 更新K值（基类通用方法）
// 更新K值（带保护）
auto Flash::updateKValues(const std::vector<double> &phi_liquid,
                          const std::vector<double> &phi_vapor)
    -> std::vector<double> {
  std::vector<double> k_current(initial_k_.size(), 0.0);
  for (size_t i = 0; i < k_current.size(); ++i)
  {
    double pv = phi_vapor[i];
    double pl = phi_liquid[i];
    if (pv < 1e-12)
      pv = 1e-12;
    if (pl < 1e-12)
      pl = 1e-12;
    k_current[i] = pl / pv;
  }
  return k_current;
}

// 多相扩展接口预留
void Flash::proposeNewPhaseByStabilityAnalysis() {
  std::cout << "稳定性分析接口预留，未实现。\n";
}

// 通用结果显示方法
void Flash::displayResults() const {
  std::cout << "闪蒸计算结果:\n";
  std::cout << "压力: " << pressure_ << " Pa\n";
  std::cout << "温度: " << temperature_ << " K\n";
  std::cout << "气相分率: " << vapor_fraction_ << "\n";

  if (!vap_comp_frac_.empty())
  {
    std::cout << "气相摩尔分率:\n";
    for (size_t i = 0; i < vap_comp_frac_.size(); ++i)
    {
      std::cout << "  组分 " << i + 1 << ": " << vap_comp_frac_[i] << "\n";
    }
  }

  if (!liq_comp_frac_.empty())
  {
    std::cout << "液相摩尔分率:\n";
    for (size_t i = 0; i < liq_comp_frac_.size(); ++i)
    {
      std::cout << "  组分 " << i + 1 << ": " << liq_comp_frac_[i] << "\n";
    }
  }
}
// 派生类 PTFlash 的构造函数 - 简化版
PTFlash::PTFlash(double pressure, double temperature,
                 const std::vector<double> &composition,
                 property_package::PropertyPackage &property_package,
                 const std::vector<double> &initial_k)
    : Flash(composition, property_package, pressure, temperature,
            std::numeric_limits<double>::quiet_NaN(), initial_k) {
  // 设置闪蒸类型
  flash_type_ = FlashType::PT;

  // 特定验证
  validateInput();

  // 如果未提供初始K值，则计算
  if (initial_k_.empty())
  {
    initializeKWithWilson(temperature_, pressure_);
  }
}

// PTFlash 的 calculate 方法实现
void PTFlash::calculate(ConvergenceMethod method) {
  calculateVaporFraction(method);
}

// PTFlash 的 calculateVaporFraction 方法
void PTFlash::calculateVaporFraction(ConvergenceMethod method) {
  std::cout << "计算PT闪蒸的气相分率...\n";

  // 设置外部迭代参数
  const int MAX_OUTER_ITERATIONS = 100;
  const double TOL_OUTER = 1e-6;
  const double TOL_INNER = 1e-6;

  // 初始化 K 值
  std::vector<double> k_current = initial_k_;
  std::vector<double> k_previous = initial_k_;

  double v_current = initial_V_; // 使用初始气相分率估计

  switch (method)
  {
  case ConvergenceMethod::NEWTON_RAPHSON: {
    for (int outer_iter = 0; outer_iter < MAX_OUTER_ITERATIONS; ++outer_iter)
    {
      std::cout << "外部迭代 " << outer_iter + 1 << ":\n";

      auto rr_func = [&](double vapfrac) -> double {
        return calculateRachfordRice(vapfrac, k_current);
      };
      auto rr_deriv = [&](double vapfrac) -> double {
        return calculateRachfordRiceDeriv(vapfrac, k_current);
      };
      // 使用牛顿-拉夫逊求解器求解 V
      double v_new = NewtonRaphsonSolver::solve(rr_func, rr_deriv, v_current,
                                                TOL_INNER, 100);
      v_new = v_new > 1 ? 0.9999 : v_new;
      v_new = v_new < 0 ? 0.0001 : v_new;
      if ((v_new == 0.9999 || v_new == 0.0001) && outer_iter > 20)
      {
        v_new = v_new == 0.9999 ? 1 : 0;
      }
      vap_comp_frac_ = computeVapCompFractions(v_new, k_current);
      liq_comp_frac_ = computeLiqCompFractions(v_new, k_current);

      // 计算液相和气相的逸度系数
      std::vector<double> phi_liquid =
          property_package_.calculateLiquidFugacityCoefficientMixture(
              temperature_, pressure_, liq_comp_frac_);
      std::vector<double> phi_vapor =
          property_package_.calculateVaporFugacityCoefficientMixture(
              temperature_, pressure_, vap_comp_frac_);

      k_current = updateKValues(phi_liquid, phi_vapor);
      // 检查收敛性
      if (checkConvergence(k_current, k_previous, v_new, v_current, TOL_OUTER,
                           outer_iter))
      {
        return;
      }

      // 更新变量以进行下一次迭代
      k_previous = k_current;
      v_current = v_new;
    }
    // 如果达到最大迭代次数仍未收敛
    throw std::runtime_error("闪蒸计算未能收敛");
  }
  case ConvergenceMethod::HALLEY: {
    for (int outer_iter = 0; outer_iter < MAX_OUTER_ITERATIONS; ++outer_iter)
    {
      // 1. 用当前K值通过Halley迭代求解新的β
      double f = calculateRachfordRice(v_current, k_current);
      double f1 = calculateRachfordRiceDeriv(v_current, k_current);
      double f2 = calculateRachfordRiceSecondDeriv(v_current, k_current);

      double denom = (2.0 * f1 * f1 - f * f2);
      if (std::abs(denom) < 1e-12)
        throw std::runtime_error("Halley分母接近零。");

      double v_new = v_current - (2.0 * f * f1) / denom;
      v_new = std::min(std::max(v_new, 1e-6), 1 - 1e-6);

      // 2. 计算新的气液相组成
      vap_comp_frac_ = computeVapCompFractions(v_new, k_current);
      liq_comp_frac_ = computeLiqCompFractions(v_new, k_current);

      // 3. 关键缺失：更新K值（根据新相组成计算逸度系数）
      std::vector<double> phi_liquid =
          property_package_.calculateLiquidFugacityCoefficientMixture(
              temperature_, pressure_, liq_comp_frac_);
      std::vector<double> phi_vapor =
          property_package_.calculateVaporFugacityCoefficientMixture(
              temperature_, pressure_, vap_comp_frac_);
      k_current = updateKValues(phi_liquid, phi_vapor); // 必须添加这一步
    std::cout << "外部迭代 " << outer_iter + 1 << ":\n";
      // 4. 检查收敛（此时k_current已更新，比较有效）
      if (checkConvergence(k_current, k_previous, v_new, v_current, TOL_OUTER,
                           outer_iter))
      {
        return;
      }

      // 5. 更新迭代变量
      v_current = v_new;
      k_previous = k_current; // 现在k_previous存储的是上一次的K值，有效
    }
    throw std::runtime_error("Halley方法未收敛");
  }
  }
}

// PTFlash 的 calculateRachfordRice 方法
auto PTFlash::calculateRachfordRice(
    const double &vapfrac, const std::vector<double> &K_CURRENT) -> double {
  double sum = 0.0;
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    double denominator = 1.0 + (vapfrac * (K_CURRENT[i] - 1.0));
    if (denominator == 0.0)
    {
      throw std::runtime_error("Rachford-Rice方程中遇到除以零。");
    }
    sum += composition_[i] * (K_CURRENT[i] - 1.0) / denominator;
  }
  return sum;
}

// PTFlash 的 calculateRachfordRiceDeriv 方法
auto PTFlash::calculateRachfordRiceDeriv(
    const double &vapfrac, const std::vector<double> &K_CURRENT) -> double {
  double sum = 0.0;
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    double denominator = 1.0 + (vapfrac * (K_CURRENT[i] - 1.0));
    if (denominator == 0.0)
    {
      throw std::runtime_error("Rachford-Rice方程导数中遇到除以零。");
    }
    sum += -composition_[i] * (K_CURRENT[i] - 1.0) * (K_CURRENT[i] - 1.0) /
           (denominator * denominator);
  }
  return sum;
}

// Rachford-Rice 二阶导数
auto PTFlash::calculateRachfordRiceSecondDeriv(
    const double &vapfrac, const std::vector<double> &K_CURRENT) -> double {
  double sum = 0.0;
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    double denominator = 1.0 + vapfrac * (K_CURRENT[i] - 1.0);
    if (denominator == 0.0)
      throw std::runtime_error("RR二阶导除零。");
    double ki_minus_1 = K_CURRENT[i] - 1.0;
    sum += 2.0 * composition_[i] * ki_minus_1 * ki_minus_1 * ki_minus_1 /
           (denominator * denominator * denominator);
  }
  return sum;
}

// PTFlash 的 checkConvergence 方法
auto PTFlash::checkConvergence(const std::vector<double> &k_current,
                               const std::vector<double> &k_previous,
                               double v_new, double v_current, double TOL_OUTER,
                               int outer_iter) -> bool {
  // 计算K值的最大变化量
  double max_k_diff = 0.0;
  for (size_t i = 0; i < k_current.size(); ++i)
  {
    double diff = std::abs(k_current[i] - k_previous[i]);
    max_k_diff = std::max(diff, max_k_diff);
  }

  // 计算气相摩尔分数的变化量
  double v_diff = std::abs(v_new - v_current);

  // 输出迭代信息，帮助监控收敛过程
  std::cout << " K的最大差异: " << max_k_diff << "\n"
            << " V的差异: " << v_diff << "\n";

  // 检查是否满足收敛条件
  if (max_k_diff < TOL_OUTER && v_diff < TOL_OUTER)
  {
    std::cout << outer_iter + 1 << " 次外部迭代后收敛。\n";

    // 更新气相分数
    vapor_fraction_ = v_new;
    return true; // 表示已收敛
  }

  return false; // 表示未收敛
}

// 派生类 PVFlash 的构造函数 - 简化版
PVFlash::PVFlash(double pressure, double vapor_fraction,
                 const std::vector<double> &composition,
                 property_package::PropertyPackage &property_package,
                 const double &initial_T, const std::vector<double> &initial_k)
    : Flash(composition, property_package, pressure,
            std::numeric_limits<double>::quiet_NaN(), vapor_fraction, initial_k,
            initial_T) {
  // 设置闪蒸类型
  flash_type_ = FlashType::PV;

  // 特定验证
  validateInput();

  // 如果未提供初始温度，则计算
  if (std::isnan(initial_T_))
  {
    initializeTWithWilson(pressure_);
    std::cout << "使用Wilson方法初始化温度: " << initial_T_ << '\n';
  }

  // 如果未提供初始K值，则计算
  if (initial_k_.empty())
  {
    initializeKWithWilson(initial_T_, pressure_);
  }
}

// PVFlash 的 calculate 方法实现
void PVFlash::calculate(ConvergenceMethod method) {
  calculateTemperature(method);
}

// PVFlash 的 calculateTemperature 方法
void PVFlash::calculateTemperature(ConvergenceMethod method) {
  std::cout << "计算PV闪蒸的温度...\n";

  // 设置迭代参数
  const int MAX_OUTER_ITERATIONS = 100;

  double T_current = initial_T_;
  double T_previous = initial_T_;
  std::vector<double> k_current = initial_k_;
  std::vector<double> k_previous = initial_k_;

  switch (method)
  {
  case ConvergenceMethod::NEWTON_RAPHSON: {
    for (int outer_iter = 0; outer_iter < MAX_OUTER_ITERATIONS; ++outer_iter)
    {
      std::cout << "迭代次数 " << outer_iter + 1 << ":\n";

      vap_comp_frac_ = computeVapCompFractions(vapor_fraction_, k_current);
      liq_comp_frac_ = computeLiqCompFractions(vapor_fraction_, k_current);

      double err =
          calculateRachfordRice(k_current, vapor_fraction_, composition_);
      std::cout << "    Rachford-Rice 方程误差: " << err << "\n";
      double derr = calculateRachfordRiceDeriv(T_current, k_current);
      std::cout << "    Rachford-Rice 方程导数: " << derr << "\n";
      T_current = updateTemperature(T_current, err, derr);
      // 使用牛顿-拉夫逊求解器求解温度 T

      std::vector<double> phi_liquid =
          property_package_.calculateLiquidFugacityCoefficientMixture(
              T_current, pressure_, liq_comp_frac_);
      std::vector<double> phi_vapor =
          property_package_.calculateVaporFugacityCoefficientMixture(
              T_current, pressure_, vap_comp_frac_);
      k_current = updateKValues(phi_liquid, phi_vapor);

      if (checkConvergence(k_current, k_previous, T_current, T_previous,
                           outer_iter))
      {
        temperature_ = T_current;
        return;
      }
      T_previous = T_current;
      k_previous = k_current;
    }
    throw std::runtime_error("PV闪蒸计算未能收敛");
  }
  case ConvergenceMethod::HALLEY: {
    const int MAX_OUTER_ITERATIONS = 100;
    double T_current = initial_T_;
    double T_previous = initial_T_;
    std::vector<double> k_current = initial_k_;
    std::vector<double> k_previous = initial_k_;

    // 便捷函数：在给定温度下评估 f(T) = RR(K(T), V, z)
    auto f_at_T = [&](double Ttest) {
      // 用当前 V、当前 K 计算一次相组成（与现有 NEWTON 分支一致的做法）
      vap_comp_frac_ = computeVapCompFractions(vapor_fraction_, k_current);
      liq_comp_frac_ = computeLiqCompFractions(vapor_fraction_, k_current);
      // 计算 φ → K(Ttest)
      auto phiL = property_package_.calculateLiquidFugacityCoefficientMixture(
          Ttest, pressure_, liq_comp_frac_);
      auto phiV = property_package_.calculateVaporFugacityCoefficientMixture(
          Ttest, pressure_, vap_comp_frac_);
      auto Ktmp = updateKValues(phiL, phiV);
      // RR 残差（内部已做 V<0.5 / V≥0.5 分段等价处理）
      return calculateRachfordRice(Ktmp, vapor_fraction_, composition_);
    };

    for (int outer_iter = 0; outer_iter < MAX_OUTER_ITERATIONS; ++outer_iter)
    {
      std::cout << "迭代次数 " << outer_iter + 1 << ":\n";

      // 1) 计算 f, f', f''（中心差分，数值稳定）
      const double eps = std::max(1e-3, std::abs(T_current) * 1e-4);
      const double f0 = f_at_T(T_current);
      const double f_p = f_at_T(T_current + eps);
      const double f_m = f_at_T(T_current - eps);
      const double f1 = (f_p - f_m) / (2.0 * eps);
      const double f2 = (f_p - 2.0 * f0 + f_m) / (eps * eps);

      // 2) Halley 步：T_new = T - 2 f f' / (2 (f')^2 - f f'')
      double denom = 2.0 * f1 * f1 - f0 * f2;
      if (std::abs(denom) < 1e-14)
      {
        // 退化时退回牛顿步
        denom = 2.0 * f1 * f1; // 相当于一步牛顿
      }
      double step = (2.0 * f0 * f1) / denom;
      double alpha = 1.0;
      double T_trial = T_current - alpha * step;

      // 3) 简单 Armijo 回溯，保证 |f| 下降（与你现有风格一致）
      double f_trial = f_at_T(T_trial);
      int backtrack = 0;
      while (std::abs(f_trial) > 0.9 * std::abs(f0) && backtrack < 8)
      {
        alpha *= 0.5;
        T_trial = T_current - alpha * step;
        f_trial = f_at_T(T_trial);
        ++backtrack;
      }

      // 下界保护，避免错误步把温度拉得过低
      T_trial = std::max(T_trial, 150.0);

      // 4) 接受新温度，刷新 φ、K、并做你的收敛判据
      T_previous = T_current;
      T_current = T_trial;

      vap_comp_frac_ = computeVapCompFractions(vapor_fraction_, k_current);
      liq_comp_frac_ = computeLiqCompFractions(vapor_fraction_, k_current);
      auto phiL = property_package_.calculateLiquidFugacityCoefficientMixture(
          T_current, pressure_, liq_comp_frac_);
      auto phiV = property_package_.calculateVaporFugacityCoefficientMixture(
          T_current, pressure_, vap_comp_frac_);
      k_previous = k_current;
      k_current = updateKValues(phiL, phiV);

      // 5) 复用你现有的收敛检查（K_maxdiff 与 ΔT）
      if (checkConvergence(k_current, k_previous, T_current, T_previous,
                           outer_iter))
      {
        temperature_ = T_current;
        return;
      }

      std::cout << "T_old : " << T_previous << "    T_new : " << T_current
                << "\n";
    }

    throw std::runtime_error("Halley方法未收敛");
  }
  }
}

// PVFlash 的 calculateRachfordRice 方法
auto PVFlash::calculateRachfordRice(
    const std::vector<double> &K_value, const double &vaporfraction,
    const std::vector<double> &composition) -> double {
  double sum = 0.0;
  double denominator = 0;
  if (vaporfraction < 0.5)
  {
    for (size_t index = 0; index < composition.size(); ++index)
    {
      denominator = 1.0 + (vaporfraction * (K_value[index] - 1.0));
      if (denominator == 0.0)
      {
        throw std::runtime_error("Rachford-Rice方程中遇到除以零。");
      }

      sum += composition[index] * (K_value[index] - 1.0) / denominator;
    }
  }
  else
  {
    for (size_t index = 0; index < composition.size(); ++index)
    {
      denominator = 1.0 + (vaporfraction * (K_value[index] - 1.0));
      if (denominator == 0.0)
      {
        throw std::runtime_error("Rachford-Rice方程中遇到除以零。");
      }
      sum += composition[index] / denominator;
    }
    sum -= 1; // 数学形式上与原rr方程是等效的，在V>0.5的情况下更稳定，误差更小
  }
  return sum;
}

// PVFlash 的 calculateRachfordRiceDeriv 方法
auto PVFlash::calculateRachfordRiceDeriv(
    const double &temperature, const std::vector<double> &K_value) -> double {
  const double V = vapor_fraction_;
  // 1) 选择稳定的差分间隔（避免 T 很小时 epsilon 过小）
  const double epsilon = std::max(1e-3, std::abs(temperature) * 1e-4);
  const double T1 = temperature - epsilon;
  const double T2 = temperature + epsilon;

  // 2) 正确的 ∂f/∂K_i（与 PVFlash::calculateRachfordRice 的分段保持一致）
  std::vector<double> dRRdKi(composition_.size(), 0.0);
  for (size_t i = 0; i < composition_.size(); ++i)
  {
    const double denom = 1.0 + V * (K_value[i] - 1.0);
    if (denom == 0.0)
      throw std::runtime_error("RR二阶/导数除零。");
    if (V < 0.5)
    {
      // f = Σ z_i (K_i-1)/(1+V(K_i-1))  ⇒ ∂f/∂K_i = z_i / denom^2
      dRRdKi[i] = composition_[i] / (denom * denom);
    }
    else
    {
      // f = Σ z_i /(1+V(K_i-1)) - 1     ⇒ ∂f/∂K_i = - z_i * V / denom^2
      dRRdKi[i] = -composition_[i] * V / (denom * denom);
    }
  }

  // 3) 差分计算 dK_i/dT（沿用你原有的 fugacity→K 流程）
  //    注意：用上一步固定的 x,y 计算 φ，即可得到 K(T±ε)
  std::vector<double> phiL1 =
      property_package_.calculateLiquidFugacityCoefficientMixture(
          T1, pressure_, liq_comp_frac_);
  std::vector<double> phiV1 =
      property_package_.calculateVaporFugacityCoefficientMixture(
          T1, pressure_, vap_comp_frac_);
  std::vector<double> phiL2 =
      property_package_.calculateLiquidFugacityCoefficientMixture(
          T2, pressure_, liq_comp_frac_);
  std::vector<double> phiV2 =
      property_package_.calculateVaporFugacityCoefficientMixture(
          T2, pressure_, vap_comp_frac_);

  std::vector<double> K1 = updateKValues(phiL1, phiV1);
  std::vector<double> K2 = updateKValues(phiL2, phiV2);

  double dRRdT = 0.0;
  for (size_t i = 0; i < K1.size(); ++i)
  {
    const double dKi_dT = (K2[i] - K1[i]) / (2.0 * epsilon);
    dRRdT += dRRdKi[i] * dKi_dT;
  }
  return dRRdT;
}

// PVFlash 的 updateTemperature 方法
auto PVFlash::updateTemperature(const double &T, const double &err,
                                const double &derr) -> double {
  if (std::abs(derr) < 1e-12)
    throw std::runtime_error("牛顿法迭代过程中分母接近0");
  double alpha = 1.0;
  double step = err / derr;

  // 原有相对步长限制先保留
  if (std::abs(step) > 0.1 * std::max(100.0, std::abs(T)))
  {
    alpha = 0.1 * std::max(100.0, std::abs(T)) / std::abs(step);
  }

  // Armijo 型回溯（需要一个计算 f 的 lambda；这里用最近一次的 k_current
  // 近似即可）
  auto f_at = [&](double Ttest) {
    // 用当前 V,k 组合出的 x,y 计算 φ → K(Ttest) → f
    auto phiL = property_package_.calculateLiquidFugacityCoefficientMixture(
        Ttest, pressure_, liq_comp_frac_);
    auto phiV = property_package_.calculateVaporFugacityCoefficientMixture(
        Ttest, pressure_, vap_comp_frac_);
    auto Ktmp = updateKValues(phiL, phiV);
    return calculateRachfordRice(Ktmp, vapor_fraction_, composition_);
  };

  const double f0 = err;
  double T_new = T - alpha * step;
  double f1 = f_at(T_new);

  int backtrack = 0;
  while (std::abs(f1) > 0.9 * std::abs(f0) && backtrack < 8)
  { // 0.9: 需要有明显下降
    alpha *= 0.5;
    T_new = T - alpha * step;
    f1 = f_at(T_new);
    ++backtrack;
  }

  // 物理下界保护，避免被错误步拉到极低温
  T_new = std::max(T_new, 150.0); // 也可用系统最小可用温度

  std::cout << "T_old : " << T << "    T_new : " << T_new << "\n";
  return T_new;
}

// PVFlash 的 checkConvergence 方法
auto PVFlash::checkConvergence(const std::vector<double> &k_current,
                               const std::vector<double> &k_previous,
                               const double &T_current,
                               const double &T_previous,
                               int outer_iter) -> bool {
  // 计算 K 值的最大变化
  double max_k_diff = 0.0;
  for (size_t i = 0; i < k_current.size(); ++i)
  {
    double diff = std::abs(k_current[i] - k_previous[i]);
    max_k_diff = std::max(max_k_diff, diff);
  }
  double diff_T = std::abs(T_current - T_previous);
  // 输出迭代信息，监控收敛过程
  std::cout << " K的最大差异: " << max_k_diff << "\n";
  std::cout << " T的差异: " << diff_T << "\n";

  // 判断是否满足收敛条件
  if (max_k_diff < 1e-6 && diff_T < 1e-3)
  {
    std::cout << outer_iter + 1 << " 次外部迭代后收敛。\n";
    return true;
  }
  return false;
}

// 派生类 TVFlash 的构造函数
TVFlash::TVFlash(double temperature, double vapor_fraction,
                 const std::vector<double> &composition,
                 property_package::PropertyPackage &property_package,
                 const double &initial_P, const std::vector<double> &initial_k)
    : Flash(composition, property_package,
            std::numeric_limits<double>::quiet_NaN(), temperature,
            vapor_fraction, initial_k, std::numeric_limits<double>::quiet_NaN(),
            initial_P) {
  // 设置闪蒸类型
  flash_type_ = FlashType::TV;

  // 特定验证
  validateInput();

  // 如果未提供初始压力，则计算
  if (std::isnan(initial_P_))
  {
    initializePWithWilson(temperature_);
    std::cout << "使用Wilson方法初始化压力: " << initial_P_ << '\n';
  }

  // 如果未提供初始K值，则计算
  if (initial_k_.empty())
  {
    initializeKWithWilson(temperature_, initial_P_);
  }
}

// TVFlash 的 calculate 方法实现
void TVFlash::calculate(ConvergenceMethod method) { calculatePressure(method); }

// TVFlash 的 calculatePressure 方法
void TVFlash::calculatePressure(ConvergenceMethod method) {
  std::cout << "计算TV闪蒸的压力...\n";

  // 设置迭代参数
  const int MAX_OUTER_ITERATIONS = 100;

  double P_current = initial_P_;
  double P_previous = initial_P_;
  std::vector<double> k_current = initial_k_;
  std::vector<double> k_previous = initial_k_;

  switch (method)
  {
  case ConvergenceMethod::NEWTON_RAPHSON: {
    for (int outer_iter = 0; outer_iter < MAX_OUTER_ITERATIONS; ++outer_iter)
    {
      std::cout << "外部迭代 " << outer_iter + 1 << ":\n";
      vap_comp_frac_ = computeVapCompFractions(vapor_fraction_, k_current);
      liq_comp_frac_ = computeLiqCompFractions(vapor_fraction_, k_current);

      double err =
          calculateRachfordRice(k_current, vapor_fraction_, composition_);
      double derr = calculateRachfordRiceDeriv(P_current, k_current);
      P_current = updatePressure(P_current, err, derr);

      std::vector<double> phi_liquid =
          property_package_.calculateLiquidFugacityCoefficientMixture(
              temperature_, P_current, liq_comp_frac_);
      std::vector<double> phi_vapor =
          property_package_.calculateVaporFugacityCoefficientMixture(
              temperature_, P_current, vap_comp_frac_);
      k_current = updateKValues(phi_liquid, phi_vapor);

      if (checkConvergence(k_current, k_previous, P_current, P_previous,
                           outer_iter))
      {
        pressure_ = P_current;
        return;
      }
      P_previous = P_current;
      k_previous = k_current;
    }
    throw std::runtime_error("TV闪蒸计算未能收敛");
  }
  case ConvergenceMethod::HALLEY: {
    auto f_at_P = [&](double Ptest){
      vap_comp_frac_ = computeVapCompFractions(vapor_fraction_, k_current);
      liq_comp_frac_ = computeLiqCompFractions(vapor_fraction_, k_current);
      auto phiL = property_package_.calculateLiquidFugacityCoefficientMixture(
          temperature_, Ptest, liq_comp_frac_);
      auto phiV = property_package_.calculateVaporFugacityCoefficientMixture(
          temperature_, Ptest, vap_comp_frac_);
      auto Ktmp = updateKValues(phiL, phiV);
      return calculateRachfordRice(Ktmp, vapor_fraction_, composition_);
    };
  
    for (int outer_iter = 0; outer_iter < MAX_OUTER_ITERATIONS; ++outer_iter) {
      std::cout << "外部迭代 " << outer_iter + 1 << ":\n";
      const double eps = std::max(1e3, std::abs(P_current) * 1e-4);
      const double f0 = f_at_P(P_current);
      const double f_p = f_at_P(P_current + eps);
      const double f_m = f_at_P(P_current - eps);
      const double f1 = (f_p - f_m) / (2.0 * eps);
      const double f2 = (f_p - 2.0*f0 + f_m) / (eps * eps);
  
      double denom = 2.0 * f1 * f1 - f0 * f2;
      if (std::abs(denom) < 1e-14) denom = 2.0 * f1 * f1; // 退回牛顿
      double step = (2.0 * f0 * f1) / denom;
      double alpha = 1.0;
      double P_trial = P_current - alpha * step;
  
      // Armijo 回溯
      double f_trial = f_at_P(P_trial);
      int bt = 0;
      while (std::abs(f_trial) > 0.9 * std::abs(f0) && bt < 8) {
        alpha *= 0.5;
        P_trial = P_current - alpha * step;
        f_trial = f_at_P(P_trial);
        ++bt;
      }
      if (P_trial <= 1.0) P_trial = std::max(0.5 * P_current, 1.0);
  
      // 刷新 φ,K 并做你的收敛判据
      P_previous = P_current;
      P_current  = P_trial;
      vap_comp_frac_ = computeVapCompFractions(vapor_fraction_, k_current);
      liq_comp_frac_ = computeLiqCompFractions(vapor_fraction_, k_current);
      auto phiL = property_package_.calculateLiquidFugacityCoefficientMixture(
          temperature_, P_current, liq_comp_frac_);
      auto phiV = property_package_.calculateVaporFugacityCoefficientMixture(
          temperature_, P_current, vap_comp_frac_);
      k_previous = k_current;
      k_current  = updateKValues(phiL, phiV);
  
      if (checkConvergence(k_current, k_previous, P_current, P_previous, outer_iter)) {
        pressure_ = P_current;
        return;
      }
    }
    throw std::runtime_error("Halley方法未收敛");
  }
  }
}

// TVFlash 的 calculateRachfordRice 方法
auto TVFlash::calculateRachfordRice(
    const std::vector<double> &K_value, const double &vaporfraction,
    const std::vector<double> &composition) -> double {
  double sum = 0.0;
  double denominator = 0;
  if (vaporfraction < 0.5)
  {
    for (size_t index = 0; index < composition.size(); ++index)
    {
      denominator = 1.0 + (vaporfraction * (K_value[index] - 1.0));
      if (denominator == 0.0)
      {
        throw std::runtime_error("Rachford-Rice方程中遇到除以零。");
      }

      sum += composition[index] * (K_value[index] - 1.0) / denominator;
    }
  }
  else
  {
    for (size_t index = 0; index < composition.size(); ++index)
    {
      denominator = 1.0 + (vaporfraction * (K_value[index] - 1.0));
      if (denominator == 0.0)
      {
        throw std::runtime_error("Rachford-Rice方程中遇到除以零。");
      }
      sum += composition[index] / denominator;
    }
    sum -= 1; // 数学形式上与原rr方程是等效的，在V>0.5的情况下更稳定，误差更小
  }
  return sum;
}

// TVFlash 的 calculateRachfordRiceDeriv 方法
auto TVFlash::calculateRachfordRiceDeriv(
  const double &pressure, const std::vector<double> &K_value) -> double {
// 数值差分步长（相对量级），并做下界保护，避免 P-ε <= 0
const double eps_rel = std::max(1e-4, 1e-6);                  // 相对步
const double epsilon = std::max(1e3, std::abs(pressure) * eps_rel); // 绝对步(>=1000 Pa)
const double P1 = std::max(1.0, pressure - epsilon);
const double P2 = pressure + epsilon;

// 计算 ∂f/∂K_i（与 RR 的分段形式严格一致）
std::vector<double> dRRdKi(composition_.size(), 0.0);
for (size_t i = 0; i < composition_.size(); ++i) {
  const double denom = 1.0 + vapor_fraction_ * (K_value[i] - 1.0);
  if (denom == 0.0) {
    throw std::runtime_error("Rachford–Rice 导数中遇到除以零。");
  }
  if (vapor_fraction_ < 0.5) {
    // f = Σ z_i (K_i-1) / (1 + V (K_i-1))  ⇒  ∂f/∂K_i = z_i / denom^2
    dRRdKi[i] = composition_[i] / (denom * denom);
  } else {
    // f = Σ z_i / (1 + V (K_i-1)) - 1      ⇒  ∂f/∂K_i = - z_i * V / denom^2
    dRRdKi[i] = -composition_[i] * vapor_fraction_ / (denom * denom);
  }
}

// 在 P±ε 处评估 φ → K，用中心差分近似 dK_i/dP
// 注意：此处使用“最近一次”的 x,y（liq_comp_frac_, vap_comp_frac_）作为固定点
// 与你现有外层迭代顺序一致。
auto phiL1 = property_package_.calculateLiquidFugacityCoefficientMixture(
    temperature_, P1, liq_comp_frac_);
auto phiV1 = property_package_.calculateVaporFugacityCoefficientMixture(
    temperature_, P1, vap_comp_frac_);
auto phiL2 = property_package_.calculateLiquidFugacityCoefficientMixture(
    temperature_, P2, liq_comp_frac_);
auto phiV2 = property_package_.calculateVaporFugacityCoefficientMixture(
    temperature_, P2, vap_comp_frac_);

std::vector<double> K1 = updateKValues(phiL1, phiV1);
std::vector<double> K2 = updateKValues(phiL2, phiV2);

// 链式法则：df/dP = Σ (∂f/∂K_i) * (dK_i/dP)
double dRRdP = 0.0;
for (size_t i = 0; i < K1.size(); ++i) {
  const double dKi_dP = (K2[i] - K1[i]) / (P2 - P1); // 中心差分
  dRRdP += dRRdKi[i] * dKi_dP;
}
return dRRdP;
}


// TVFlash 的 updatePressure 方法
auto TVFlash::updatePressure(const double &P, const double &err, const double &derr) -> double {
  if (std::abs(derr) < 1e-12) throw std::runtime_error("牛顿法迭代过程中分母接近0");
  double step = err / derr;
  double alpha = 1.0;

  // 相对步长硬限（保留）
  if (std::abs(step) > 0.1 * std::max(1e5, std::abs(P))) {
    alpha = 0.1 * std::max(1e5, std::abs(P)) / std::abs(step);
  }

  // 定义 f(P) 评估（用当前 x,y 计算 φ→K，再算 RR）
  auto f_at = [&](double Ptest) {
    auto phiL = property_package_.calculateLiquidFugacityCoefficientMixture(
        temperature_, Ptest, liq_comp_frac_);
    auto phiV = property_package_.calculateVaporFugacityCoefficientMixture(
        temperature_, Ptest, vap_comp_frac_);
    auto Ktmp = updateKValues(phiL, phiV);
    return calculateRachfordRice(Ktmp, vapor_fraction_, composition_);
  };

  const double f0 = err;
  double P_new = P - alpha * step;
  double f1 = f_at(P_new);
  int backtrack = 0;
  while (std::abs(f1) > 0.9 * std::abs(f0) && backtrack < 8) {
    alpha *= 0.5;
    P_new = P - alpha * step;
    f1 = f_at(P_new);
    ++backtrack;
  }
  if (P_new <= 1.0) P_new = std::max(0.5 * P, 1.0); // 物理下界保护

  std::cout << "P_old : " << P << "    P_new : " << P_new << "\n";
  return P_new;
}


// TVFlash 的 checkConvergence 方法
auto TVFlash::checkConvergence(const std::vector<double> &k_current,
                               const std::vector<double> &k_previous,
                               const double &P_current,
                               const double &P_previous,
                               int outer_iter) -> bool {
  // 计算 K 值的最大变化
  double max_k_diff = 0.0;
  for (size_t i = 0; i < k_current.size(); ++i)
  {
    double diff = std::abs(k_current[i] - k_previous[i]);
    max_k_diff = std::max(max_k_diff, diff);
  }
  double diff_P = std::abs(P_current - P_previous) / P_current;
  // 输出迭代信息，监控收敛过程
  std::cout << " K的最大差异: " << max_k_diff << "\n";
  std::cout << " P的相对差异: " << diff_P << "\n";

  // 判断是否满足收敛条件
  if (max_k_diff < 1e-6 && diff_P < 1e-3)
  {
    std::cout << outer_iter + 1 << " 次外部迭代后收敛。\n";
    return true;
  }
  return false;
}