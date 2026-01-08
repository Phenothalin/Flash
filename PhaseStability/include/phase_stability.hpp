#pragma once

#include "thermo_backend.hpp"

#include <vector>
#include <string>

namespace phase_stability {

struct StabilityOptions {
  int max_ss_iters = 100;          // Michelsen successive substitution 最大迭代次数
  int n_random_seeds = 8;          // 随机种子数量（用于捕捉多个极小点）
  double tpd_tol = 1e-12;          // 判定不稳定阈值：TPD < -tpd_tol
  double comp_floor = 1e-14;       // 组成下限（避免 ln(0)）
  double distinct_l1 = 1e-2;       // 判定两个解“不同”的 L1 距离阈值
  int max_incipient = 4;        // 最多保留的 incipient 数量（用于 LLV/LLE 初始化）
  bool dual_reference = true;   // 同时以 V-like 与 L-like 作为参考相做稳定性分析
  bool verbose = false;
};

struct IncipientPhase {
  int phase_flag = 0;              // thermo backend phaseFlag
  int reference_phase_flag = 0;    // 本次 TPD/切平面所使用的参考相 phaseFlag
  std::vector<double> x;           // trial phase composition (mole fraction)
  double tpd = 0.0;                // tangent plane distance at optimum
  int iters = 0;
};

struct StabilityResult {
  bool stable = true;              // 对参考单相解是否稳定
  int reference_phase_flag = 0;    // 稳定性测试使用的参考相标志
  double g_ref = 0.0;              // 参考相的摩尔 Gibbs 能（用于单相判别）
  std::vector<IncipientPhase> incipient; // 发现的“起始相/分裂相”候选（去重后）
};

class PhaseStabilityAnalyzer {
public:
  explicit PhaseStabilityAnalyzer(const thermo::IThermoBackend& thermo);

  // 输入 z：可为摩尔分率或摩尔数（内部会归一化）
  StabilityResult analyze(double T, double P, const std::vector<double>& z,
                          const StabilityOptions& opt = {}) const;

private:
  const thermo::IThermoBackend& thermo_;

  static std::vector<double> normalize_positive(const std::vector<double>& v, double floor);
  static double l1_dist(const std::vector<double>& a, const std::vector<double>& b);

  // 从 mu 得到 lnphi：mu = RT*(lnphi + ln x + ln P)
  static std::vector<double> lnphi_from_mu(double T, double P,
                                           const std::vector<double>& x,
                                           const std::vector<double>& mu);

  static double molar_gibbs(const std::vector<double>& x, const std::vector<double>& mu);

  std::vector<std::vector<double>> build_seeds(double T, double P,
                                               const std::vector<double>& z_norm,
                                               int trial_phase_flag,
                                               const StabilityOptions& opt) const;

  bool successive_substitution(double T, double P,
                               const std::vector<double>& z_norm,
                               const std::vector<double>& lnphi_z,
                               int trial_phase_flag,
                               std::vector<double>& w_inout,
                               int& iters,
                               const StabilityOptions& opt,
                               std::vector<double>* lnphi_w_out = nullptr) const;
};

} // namespace phase_stability
