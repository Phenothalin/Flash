#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <cmath>

static constexpr double R_CONST = 8.314462618;

// -------- ThermoPack C++ headers (兼容两种目录拼写) --------
#if __has_include(<cppThermopack/cubic.h>)
  #include <cppThermopack/cubic.h>
  #include <cppThermopack/thermo.h>
#elif __has_include(<cppThermoPack/cubic.h>)
  #include <cppThermoPack/cubic.h>
  #include <cppThermoPack/thermo.h>
#else
  #error "Cannot find ThermoPack C++ headers: cppThermopack/cubic.h"
#endif

namespace thermo {  // randflash namespace (随你项目调整)

// 小工具：把相内摩尔数 n 归一为 x，并返回总摩尔数 N
inline std::pair<std::vector<double>, double>
normalize_n(const std::vector<double>& n) {
  double N = std::accumulate(n.begin(), n.end(), 0.0);
  if (N <= 0.0) throw std::runtime_error("normalize_n: total moles <= 0");
  std::vector<double> x(n.size());
  for (size_t i = 0; i < n.size(); ++i) x[i] = n[i] / N;
  return {x, N};
}

// ThermoPack 适配器（立方型 EOS 为例；构造时选择 PR/SRK 等）
class ThermoAdapterTP {
public:
  // 构造：components 为 ThermoPack 识别的组分字符串（逗号分隔）
  // eos: "PR"/"SRK"；mixing: "vdW"/...；alpha: "Classic"/...；ref: "Default"/...
  explicit ThermoAdapterTP(const std::string& components_csv,
                           const std::string& eos   = "PR",
                           const std::string& mixing= "vdW",
                           const std::string& alpha = "Classic",
                           const std::string& ref   = "Default",
                           bool volume_shift = false)
  : eos_(components_csv, eos, mixing, alpha, ref, volume_shift),
    liqph_(eos_.LIQPH),
    vapph_(eos_.VAPPH),
    mingibbsph_(4) {}  // Phase::mingibbs = 4 in ThermoPack

  // 基本信息
  int LIQPH() const { return liqph_; }
  int VAPPH() const { return vapph_; }
  int MINGIBBSPH() const { return mingibbsph_; }  // 自动选择Gibbs能最低的根

  // -------- 1) d(ln phi)/dn：来自 ThermoPack 的公开导数 --------
  // x = n/sum(n)；调用 TP: thermo(T,P,x,phase, ..., dlnfugdn=true).dn()
  std::vector<std::vector<double>>
  dlnphi_dn(double T, double P,
            const std::vector<double>& x_or_n, // 接受 x 或 n：会自动识别
            int phase) const
  {
    std::vector<double> x = x_or_n;
    double total_moles = 1.0;
    // 简单判断：若和为 1±1e-12 视为 x；否则按 n 归一
    double s = std::accumulate(x.begin(), x.end(), 0.0);
    if (std::abs(s - 1.0) > 1e-8) {
      auto xN = normalize_n(x_or_n);
      x = std::move(xN.first);
      total_moles = xN.second;
    }

    // ThermoPack's derivative here is evaluated on the normalized composition.
    // Convert it to a derivative with respect to phase mole numbers n by
    // applying the chain-rule scale factor dx/dn ~ 1/N.
    auto prop = eos_.thermo(T, P, x, phase, /*dlnfugdt*/false,
                                      /*dlnfugdp*/false,
                                      /*dlnfugdn*/true);
    auto dlnphi = prop.dn(); // nc x nc
    if (total_moles > 0.0 && std::abs(total_moles - 1.0) > 1e-12) {
      const double inv_total = 1.0 / total_moles;
      for (auto& row : dlnphi) {
        for (double& v : row) {
          v *= inv_total;
        }
      }
    }
    return dlnphi;
  }

  // -------- 2) 组装 d(ln f)/dn = d(ln phi)/dn + d(ln x)/dn --------
  // 注意：ln p 对 n 的导数为 0（相内摩尔数变化不改 p）
  std::vector<std::vector<double>>
  assemble_dlnf_dn(double T, double P,
                   const std::vector<double>& n,
                   int phase) const
  {
    const auto& xN = normalize_n(n);
    const std::vector<double>& x = xN.first;
    const double N = xN.second;

    // dlnphi_dn must be evaluated against the phase mole numbers n, not the
    // normalized composition x. Otherwise the 1/N chain-rule scaling is lost,
    // while d(ln x)/dn below is still built in mole-number coordinates.
    auto dlnphi = dlnphi_dn(T, P, n, phase);
    const size_t nc = x.size();

    // d(ln x_i)/dn_j = (δ_ij / n_i) - 1/N
    std::vector<std::vector<double>> J(nc, std::vector<double>(nc, 0.0));
    for (size_t i = 0; i < nc; ++i) {
      const double ni = std::max(n[i], 1e-300);
      for (size_t j = 0; j < nc; ++j) {
        const double dlnx = ((i == j) ? 1.0/ni : 0.0) - 1.0/N;
        J[i][j] = dlnphi[i][j] + dlnx;
      }
    }
    return J;
  }

  // -------- 3) d(mu)/dn = RT * d(ln f)/dn --------
  std::vector<std::vector<double>>
  dmu_dn(double T, double P,
         const std::vector<double>& n,
         int phase) const
  {
    const double RT = R_CONST * T;
    auto J = assemble_dlnf_dn(T, P, n, phase);
    for (auto& row : J) for (auto& v : row) v *= RT;
    return J;
  }

  // -------- 4) 化学势 μ：使用 thermopack 的 chemical_potential_tv (理想+剩余) --------
  // 对非反应相平衡求解，modified RAND 实际需要的是与 ln(f_i) 等价的势。
  // 直接由 ln(phi_i) + ln(x_i) + ln(P) 组装可避免对总相摩尔数 N 的伪依赖。
  std::vector<double> chemicalPotentials(double T, double P,
             const std::vector<double>& n_or_x,
             int phase) const
  {
    auto nx = normalize_n(n_or_x);
    const auto& x = nx.first;
    auto prop = eos_.thermo(T, P, x, phase,
                            /*dlnfugdt*/false,
                            /*dlnfugdp*/false,
                            /*dlnfugdn*/false);
    const auto& lnphi = prop.value();

    const double RT = R_CONST * T;
    std::vector<double> mu(x.size(), 0.0);
    for (size_t i = 0; i < x.size(); ++i) {
      mu[i] = RT * (lnphi[i] + std::log(std::max(x[i], 1e-300)) + std::log(P));
    }
    return mu;
  }
  
 // -------- 5) Wilson K 值 --------
  void wilsonK(double T, double P, std::vector<double>& K) const{
    eos_.wilsonK(T, P, K);
  }

  // -------- 6) 逸度系数 φ：从化学势计算 --------
  // ln(φ_i) = μ_i/(RT) - ln(x_i) - ln(P)
  std::vector<double> lnFugacityCoefficients(double T, double P,
                                              const std::vector<double>& n_or_x,
                                              int phase) const
  {
    auto nx = normalize_n(n_or_x);
    const auto& x = nx.first;
    const size_t nc = x.size();

    // 获取化学势
    auto mu = chemicalPotentials(T, P, n_or_x, phase);

    const double RT = R_CONST * T;
    std::vector<double> lnphi(nc);

    for (size_t i = 0; i < nc; ++i) {
      // ln(φ_i) = μ_i/(RT) - ln(x_i) - ln(P)
      // 注意：P的单位需要与化学势一致，ThermoPack使用Pa
      lnphi[i] = mu[i] / RT - std::log(std::max(x[i], 1e-300)) - std::log(P);
    }

    return lnphi;
  }

  // -------- 7) 逸度系数 φ（非对数形式）--------
  std::vector<double> fugacityCoefficients(double T, double P,
                                            const std::vector<double>& n_or_x,
                                            int phase) const
  {
    auto lnphi = lnFugacityCoefficients(T, P, n_or_x, phase);
    std::vector<double> phi(lnphi.size());
    for (size_t i = 0; i < lnphi.size(); ++i) {
      phi[i] = std::exp(lnphi[i]);
    }
    return phi;
  }

  // -------- 8) 压缩因子 Z = PVm/(RT) --------
  double compressibilityFactor(double T, double P,
                               const std::vector<double>& n_or_x,
                               int phase) const
  {
    auto nx = normalize_n(n_or_x);
    const auto& x = nx.first;

    // 获取摩尔体积 Vm (m³/mol)
    auto v_prop = eos_.specific_volume(T, P, x, phase);
    double Vm = v_prop.value();

    // Z = PVm/(RT)
    return P * Vm / (R_CONST * T);
  }

private:
  // 直接用 ThermoPack 的 Cubic 封装（构造时已选定 PR/SRK 等）
  Cubic eos_;
  int liqph_{1}, vapph_{2}, mingibbsph_{4};
};

} // namespace rf
