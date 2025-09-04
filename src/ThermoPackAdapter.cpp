#include "thermo/ThermoPackAdapter.hpp"
#include <cmath>
#include <algorithm>

// ThermoPack 头文件（根据实际安装位置，通常为如下命名空间/路径）
#include "cppThermoPack/thermo.h"
#include "cppThermoPack/cubic.h"
// 如果你选 PCSAFT，则包含 pcsaft.h；下同

namespace thermo {

struct ThermoPackAdapter::Impl {
  std::unique_ptr<tp::Thermo> eos;  // 基类指针

  Impl(const std::string& eos_name, const std::vector<std::string>& comps) {
    // 例：只示范 cubic；如需 PCSAFT 改为 tp::PCSAFT(...)
    if (eos_name == "SRK" || eos_name == "PR") {
      eos = std::make_unique<tp::Cubic>(eos_name, comps);
    } else {
      // 也可扩展到 PCSAFT/SAFT-VR-Mie
      throw std::invalid_argument("Unsupported EOS in demo: " + eos_name);
    }
  }
};

ThermoPackAdapter::ThermoPackAdapter(const std::string& eos,
                                     const std::vector<std::string>& comps)
  : impl_(std::make_unique<Impl>(eos, comps))
{}

std::vector<double> ThermoPackAdapter::chemicalPotentials(const PhaseState& s) {
  // ThermoPack 接口：给定 T,P,x 计算 μ_i
  // 这里 s.n 是该相的摩尔数，先转 x
  double beta = 0.0;
  for (double v : s.n) beta += v;
  if (beta <= 0.0) throw std::invalid_argument("Phase mole total <= 0");
  std::vector<double> x(s.n.size());
  for (size_t i=0;i<x.size();++i) x[i] = s.n[i]/beta;

  // 取化学势（J/mol）
  // 假设 tp::Thermo 暴露 chemical_potential(T,P,x)；若名称不同，按实际 API 替换
  return impl_->eos->chemical_potential(s.T, s.P, x);
}

std::vector<std::vector<double>>
ThermoPackAdapter::dMu_dN(const PhaseState& s) {
  // 数值差分：对每个 n_k 做微小扰动，计算 μ_i 的变化
  const size_t C = s.n.size();
  std::vector<std::vector<double>> J(C, std::vector<double>(C, 0.0));

  // 基点 μ
  auto mu0 = chemicalPotentials(s);

  for (size_t k=0;k<C;++k) {
    PhaseState sp = s;
    double nk = std::max(1e-12, sp.n[k]);
    double h  = std::max(1e-12, std::abs(nk)*rel_eps_);
    sp.n[k] += h;

    auto mu1 = chemicalPotentials(sp);
    for (size_t i=0;i<C;++i) {
      J[i][k] = (mu1[i] - mu0[i]) / h;  // [J/mol^2]
    }
  }
  return J;
}

} // namespace thermo
