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
    vapph_(eos_.VAPPH) {}

  // 基本信息
  int LIQPH() const { return liqph_; }
  int VAPPH() const { return vapph_; }

  // -------- 1) d(ln phi)/dn：来自 ThermoPack 的公开导数 --------
  // x = n/sum(n)；调用 TP: thermo(T,P,x,phase, ..., dlnfugdn=true).dn()
  std::vector<std::vector<double>>
  dlnphi_dn(double T, double P,
            const std::vector<double>& x_or_n, // 接受 x 或 n：会自动识别
            int phase) const
  {
    std::vector<double> x = x_or_n;
    // 简单判断：若和为 1±1e-12 视为 x；否则按 n 归一
    double s = std::accumulate(x.begin(), x.end(), 0.0);
    if (std::abs(s - 1.0) > 1e-8) {
      x = normalize_n(x_or_n).first;
    }
    /* dlnfugdn在∑n = const这一约束下数值等价于dlnfugdx(易推导)
    因此后续构建的dmudn、m矩阵都是在∑n = const的约束下的（等价于gibbs-duhem约束）*/
    auto prop = eos_.thermo(T, P, x, phase, /*dlnfugdt*/false,
                                      /*dlnfugdp*/false,
                                      /*dlnfugdn*/true);
    return prop.dn(); // nc x nc
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

    auto dlnphi = dlnphi_dn(T, P, x, phase);
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

  // -------- 4) ln phi & mu/RT：仅在你确认有“值的公开 getter”时启用 --------
  std::vector<double>chemicalPotentials(double T, double P,
             const std::vector<double>& n_or_x,
             int phase) const
  {
    // std::cout<<"调用chemicalPotentials前 phaseFlag"<<phase<<std::endl;
    auto nx = normalize_n(n_or_x);
    const auto& x = nx.first;
    const double RT = R_CONST * T;
    // 1) 拿 lnφ 的“值向量”
    auto prop = eos_.thermo(T, P, x, phase, /*dlnfugdt*/false,
                                      /*dlnfugdp*/false,
                                      /*dlnfugdn*/false);
    const std::vector<double>& lnphi = prop.value();
    // for(size_t i = 0; i < lnphi.size() ; ++i){
    //   std::cout<<lnphi[i]<<"  ";
    // }
    // NOTE: avoid spamming stdout in tight loops
    // 2) μ = RT * ( ln f = lnφ + ln x + ln p )
    std::vector<double> muRT(lnphi.size());
    const double lnP = std::log(P);
    for (size_t i=0;i<muRT.size();++i)
    muRT[i] = lnphi[i] + std::log(std::max(x[i],1e-300)) + lnP;

    std::vector<double> mu(muRT.size());
    for(size_t i=0;i<muRT.size();++i)
      mu[i] = RT * muRT[i];
    return mu;
  }
  
 // -------- 5) Wilson K 值 --------
  void wilsonK(double T, double P, std::vector<double>& K) const{
    eos_.wilsonK(T, P, K);
  }
private:
  // 直接用 ThermoPack 的 Cubic 封装（构造时已选定 PR/SRK 等）
  Cubic eos_;
  int liqph_{1}, vapph_{2};
};

} // namespace rf
