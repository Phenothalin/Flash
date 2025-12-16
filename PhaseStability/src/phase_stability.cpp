#include "phase_stability.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>
#include <iostream>

namespace phase_stability {

static constexpr double R_CONST = 8.314462618;

PhaseStabilityAnalyzer::PhaseStabilityAnalyzer(const thermo::IThermoBackend& thermo)
  : thermo_(thermo) {}

std::vector<double> PhaseStabilityAnalyzer::normalize_positive(const std::vector<double>& v, double floor)
{
  if (v.empty()) return {};
  std::vector<double> x = v;
  for (double& xi : x) {
    if (!std::isfinite(xi) || xi < 0.0) xi = 0.0;
    if (xi < floor) xi = floor;
  }
  double s = std::accumulate(x.begin(), x.end(), 0.0);
  if (s <= 0.0) {
    const double uni = 1.0 / static_cast<double>(x.size());
    std::fill(x.begin(), x.end(), uni);
    return x;
  }
  for (double& xi : x) xi /= s;
  return x;
}

double PhaseStabilityAnalyzer::l1_dist(const std::vector<double>& a, const std::vector<double>& b)
{
  if (a.size() != b.size()) return 1e100;
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) s += std::abs(a[i] - b[i]);
  return s;
}

std::vector<double> PhaseStabilityAnalyzer::lnphi_from_mu(double T, double P,
                                                          const std::vector<double>& x,
                                                          const std::vector<double>& mu)
{
  const double RT = R_CONST * T;
  const double lnP = std::log(std::max(P, 1e-300));
  if (mu.size() != x.size()) {
    throw std::runtime_error("lnphi_from_mu: size mismatch");
  }
  std::vector<double> lnphi(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    const double xi = std::max(x[i], 1e-300);
    lnphi[i] = (mu[i] / RT) - std::log(xi) - lnP;
  }
  return lnphi;
}

double PhaseStabilityAnalyzer::molar_gibbs(const std::vector<double>& x, const std::vector<double>& mu)
{
  if (x.size() != mu.size()) throw std::runtime_error("molar_gibbs: size mismatch");
  double g = 0.0;
  for (size_t i = 0; i < x.size(); ++i) g += x[i] * mu[i];
  return g;
}

std::vector<std::vector<double>> PhaseStabilityAnalyzer::build_seeds(
    double T, double P,
    const std::vector<double>& z_norm,
    int trial_phase_flag,
    const StabilityOptions& opt) const
{
  const size_t C = z_norm.size();
  std::vector<std::vector<double>> seeds;
  if (C == 0) return seeds;

  // seed 0: z itself
  seeds.push_back(z_norm);

  // Wilson K-based seeds
  std::vector<double> K(C, 1.0);
  try {
    thermo_.wilsonK(T, P, K);
  } catch (...) {
    // ignore
  }

  std::vector<double> zk(C), z_over_k(C);
  for (size_t i = 0; i < C; ++i) {
    const double Ki = std::max(K[i], 1e-50);
    zk[i] = z_norm[i] * Ki;
    z_over_k[i] = z_norm[i] / Ki;
  }
  if (trial_phase_flag == thermo_.vaporPhaseFlag()) {
    seeds.push_back(normalize_positive(zk, opt.comp_floor));
  } else {
    seeds.push_back(normalize_positive(z_over_k, opt.comp_floor));
  }

  // Extreme component-emphasis seeds (helpful for LLE / immiscible behavior)
  // Pick up to 2 dominant components
  std::vector<size_t> idx(C);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return z_norm[a] > z_norm[b]; });
  const int n_dom = std::min<int>(2, static_cast<int>(C));
  for (int k = 0; k < n_dom; ++k) {
    size_t imax = idx[k];
    std::vector<double> w = z_norm;
    // put 90% into one component, remaining 10% proportional to z
    const double rem = 0.10;
    for (size_t i = 0; i < C; ++i) w[i] = rem * z_norm[i];
    w[imax] += 0.90;
    seeds.push_back(normalize_positive(w, opt.comp_floor));
  }

  // If water exists, add a nearly-pure-water seed (capture aqueous phase)
  auto names = thermo_.getComponentNames();
  for (size_t i = 0; i < std::min(names.size(), C); ++i) {
    std::string n = names[i];
    std::transform(n.begin(), n.end(), n.begin(), ::toupper);
    if (n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
      std::vector<double> w(C, opt.comp_floor);
      w[i] = 1.0;
      seeds.push_back(normalize_positive(w, opt.comp_floor));
      break;
    }
  }

  // Random Dirichlet-like seeds (biased by z)
  std::mt19937_64 rng(1234567);
  std::gamma_distribution<double> gamma(1.0, 1.0);
  for (int r = 0; r < opt.n_random_seeds; ++r) {
    std::vector<double> w(C, 0.0);
    for (size_t i = 0; i < C; ++i) {
      const double g = gamma(rng);
      // bias towards z: larger z gets larger expected w
      w[i] = g * (0.2 + 0.8 * z_norm[i]);
    }
    seeds.push_back(normalize_positive(w, opt.comp_floor));
  }

  return seeds;
}

bool PhaseStabilityAnalyzer::successive_substitution(
    double T, double P,
    const std::vector<double>& z_norm,
    const std::vector<double>& lnphi_z,
    int trial_phase_flag,
    std::vector<double>& w_inout,
    int& iters,
    const StabilityOptions& opt,
    std::vector<double>* lnphi_w_out) const
{
  const size_t C = z_norm.size();
  const double lnP = std::log(std::max(P, 1e-300));

  std::vector<double> w = normalize_positive(w_inout, opt.comp_floor);
  std::vector<double> lnphi_w(C, 0.0);

  for (iters = 0; iters < opt.max_ss_iters; ++iters) {
    // mu = RT*(lnphi + ln w + lnP)
    thermo::PhaseState st{T, P, w, trial_phase_flag};
    const std::vector<double> mu = thermo_.chemicalPotentials(st);
    lnphi_w = lnphi_from_mu(T, P, w, mu);

    std::vector<double> w_new(C, 0.0);
    for (size_t i = 0; i < C; ++i) {
      // w_i^{new} ∝ z_i * exp( lnphi_z_i - lnphi_w_i )
      const double zi = std::max(z_norm[i], opt.comp_floor);
      const double ex = std::exp(std::clamp(lnphi_z[i] - lnphi_w[i], -80.0, 80.0));
      w_new[i] = zi * ex;
    }
    w_new = normalize_positive(w_new, opt.comp_floor);

    double diff = 0.0;
    for (size_t i = 0; i < C; ++i) diff = std::max(diff, std::abs(w_new[i] - w[i]));
    w.swap(w_new);
    if (diff < 1e-10) {
      ++iters;
      break;
    }
  }

  w_inout = w;
  if (lnphi_w_out) *lnphi_w_out = lnphi_w;
  (void)lnP;
  return true;
}

StabilityResult PhaseStabilityAnalyzer::analyze(double T, double P, const std::vector<double>& z,
                                               const StabilityOptions& opt) const
{
  StabilityResult out;
  if (z.empty()) return out;

  const std::vector<double> z_norm = normalize_positive(z, opt.comp_floor);

  // Choose reference phase: compare molar Gibbs in vapor vs liquid at same z
  const int vap = thermo_.vaporPhaseFlag();
  const int liq = thermo_.liquidPhaseFlag();

  const std::vector<double> mu_v = thermo_.chemicalPotentials(thermo::PhaseState{T, P, z_norm, vap});
  const std::vector<double> mu_l = thermo_.chemicalPotentials(thermo::PhaseState{T, P, z_norm, liq});
  const double g_v = molar_gibbs(z_norm, mu_v);
  const double g_l = molar_gibbs(z_norm, mu_l);

  int ref_flag = (g_l <= g_v) ? liq : vap;
  out.reference_phase_flag = ref_flag;
  out.g_ref = (g_l <= g_v) ? g_l : g_v;

  const std::vector<double> mu_ref = (ref_flag == liq) ? mu_l : mu_v;
  const std::vector<double> lnphi_z = lnphi_from_mu(T, P, z_norm, mu_ref);

  struct Candidate {
    int phase_flag;
    std::vector<double> x;
    double tpd;
    int iters;
  };
  std::vector<Candidate> cands;

  auto try_phase = [&](int trial_flag) {
    auto seeds = build_seeds(T, P, z_norm, trial_flag, opt);
    for (auto w0 : seeds) {
      std::vector<double> w = normalize_positive(w0, opt.comp_floor);
      int iters = 0;
      std::vector<double> lnphi_w;
      successive_substitution(T, P, z_norm, lnphi_z, trial_flag, w, iters, opt, &lnphi_w);

      // Compute TPD at w
      double tpd = 0.0;
      for (size_t i = 0; i < w.size(); ++i) {
        const double wi = std::max(w[i], opt.comp_floor);
        const double zi = std::max(z_norm[i], opt.comp_floor);
        tpd += wi * (std::log(wi) + lnphi_w[i] - std::log(zi) - lnphi_z[i]);
      }

      if (std::isfinite(tpd) && tpd < -opt.tpd_tol) {
        cands.push_back(Candidate{trial_flag, w, tpd, iters});
      }
    }
  };

  // Test both vapor-like and liquid-like trial phases against the chosen reference tangent plane
  try_phase(vap);
  try_phase(liq);

  // Deduplicate by (phase_flag, composition) L1 distance
  std::vector<Candidate> uniq;
  for (const auto& c : cands) {
    bool merged = false;
    for (auto& u : uniq) {
      if (u.phase_flag == c.phase_flag && l1_dist(u.x, c.x) < opt.distinct_l1) {
        if (c.tpd < u.tpd) u = c;
        merged = true;
        break;
      }
    }
    if (!merged) uniq.push_back(c);
  }

  std::sort(uniq.begin(), uniq.end(), [](const Candidate& a, const Candidate& b){
    return a.tpd < b.tpd;
  });

  out.stable = uniq.empty();
  out.incipient.clear();
  for (const auto& u : uniq) {
    out.incipient.push_back(IncipientPhase{u.phase_flag, u.x, u.tpd, u.iters});
    if (out.incipient.size() >= 2) break; // 最多保留两个：用于 3 相初始化
  }

  if (opt.verbose) {
    std::cout << "[PhaseStability] ref=" << (ref_flag == liq ? "LIQ" : "VAP")
              << " g_ref=" << out.g_ref
              << " stable=" << out.stable
              << " n_incipient=" << out.incipient.size() << "\n";
    for (size_t k = 0; k < out.incipient.size(); ++k) {
      std::cout << "  cand[" << k << "] flag="
                << (out.incipient[k].phase_flag == liq ? "LIQ" : "VAP")
                << " tpd=" << out.incipient[k].tpd << "\n";
    }
  }

  return out;
}

} // namespace phase_stability
