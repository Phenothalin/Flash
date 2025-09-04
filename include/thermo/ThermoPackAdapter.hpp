#pragma once
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

namespace thermo {

struct PhaseState {
  double T;           // K
  double P;           // Pa
  std::vector<double> n; // mole numbers of components in THIS phase
};

class ThermoPackAdapter {
public:
  // comps: e.g. {"methane","ethane"} ; eos: "SRK" | "PR" | "PCSAFT" ...
  ThermoPackAdapter(const std::string& eos,
                    const std::vector<std::string>& comps);

  void setPhaseIndex(int j) { phase_index_ = j; }

  // chemical potentials μ_i [J/mol]
  std::vector<double> chemicalPotentials(const PhaseState& state);

  // Jacobian dμ_i/dn_k [J/mol^2] with finite-difference fallback
  std::vector<std::vector<double>> dMu_dN(const PhaseState& state);

private:
  // pimpl 指针，内部持有 ThermoPack C++ 对象（不同 EOS 的基类指针）
  struct Impl;
  std::unique_ptr<Impl> impl_;
  int phase_index_{0};

  // 数值差分参数
  double rel_eps_ = 1e-7;
};

} // namespace thermo
