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

  const int vaporFlag  = thermo_.vaporPhaseFlag();
  const int liquidFlag = thermo_.liquidPhaseFlag();

  // Normalize z (accepts mole numbers or mole fractions)
  const std::vector<double> z_norm = normalize_positive(z, opt.comp_floor);

  // Compute V-like and L-like single-phase chemical potentials at z
  const std::vector<double> mu_v = thermo_.chemicalPotentials(thermo::PhaseState{T, P, z_norm, vaporFlag});
  const std::vector<double> mu_l = thermo_.chemicalPotentials(thermo::PhaseState{T, P, z_norm, liquidFlag});

  const double g_v = molar_gibbs(z_norm, mu_v);
  const double g_l = molar_gibbs(z_norm, mu_l);

  // For reporting only: choose the lower-Gibbs single-phase branch as reference.
  const int report_ref = (g_l <= g_v) ? liquidFlag : vaporFlag;
  out.reference_phase_flag = report_ref;
  out.g_ref = (g_l <= g_v) ? g_l : g_v;

  struct Cand {
    int phase_flag;
    int ref_flag;
    std::vector<double> x;
    double tpd;
    int iters;
  };

  std::vector<Cand> all;

  auto run_with_reference = [&](int ref_flag, const std::vector<double>& mu_ref) {
    const std::vector<double> lnphi_z_ref = lnphi_from_mu(T, P, z_norm, mu_ref);

    auto try_phase = [&](int trial_phase_flag) {
      auto seeds = build_seeds(T, P, z_norm, trial_phase_flag, opt);
      for (auto w0 : seeds) {
        std::vector<double> w = normalize_positive(w0, opt.comp_floor);
        int iters = 0;
        std::vector<double> lnphi_w;

        (void)successive_substitution(
            T, P, z_norm, lnphi_z_ref, trial_phase_flag,
            w, iters, opt, &lnphi_w);

        // TPD = Σ w_i [ ln(w_i) + lnphi(w) - ln(z_i) - lnphi_ref(z) ]
        double tpd = 0.0;
        for (size_t i = 0; i < w.size(); ++i) {
          tpd += w[i] * ((std::log(w[i]) + lnphi_w[i])
                       - (std::log(z_norm[i]) + lnphi_z_ref[i]));
        }

        if (tpd < -opt.tpd_tol) {
          all.push_back(Cand{trial_phase_flag, ref_flag, w, tpd, iters});
        }
      }
    };

    // Try both trial phase types against this reference tangent plane.
    try_phase(vaporFlag);
    try_phase(liquidFlag);
  };

  if (opt.dual_reference) {
    run_with_reference(vaporFlag, mu_v);
    run_with_reference(liquidFlag, mu_l);
  } else {
    run_with_reference(report_ref, (report_ref == liquidFlag) ? mu_l : mu_v);
  }

  // ------------------------
  //  Incipient clustering / de-dup
  // ------------------------
  std::sort(all.begin(), all.end(), [](const Cand& a, const Cand& b) {
    return a.tpd < b.tpd;
  });

  std::vector<Cand> uniq;
  for (const auto& c : all) {
    bool merged = false;
    for (auto& u : uniq) {
      if (u.phase_flag != c.phase_flag) continue;
      if (l1_dist(u.x, c.x) < opt.distinct_l1) {
        // Keep the more negative TPD (stronger instability)
        if (c.tpd < u.tpd) u = c;
        merged = true;
        break;
      }
    }
    if (!merged) uniq.push_back(c);
  }

  out.stable = uniq.empty();
  if (out.stable) return out;

  // Prefer returning: (V + L) and if possible two distinct liquids (for LLE/LLV init)
  std::vector<Cand> vap, liq;
  for (const auto& u : uniq) {
    if (u.phase_flag == vaporFlag) vap.push_back(u);
    else liq.push_back(u);
  }
  auto sort_tpd = [](const Cand& a, const Cand& b){ return a.tpd < b.tpd; };
  std::sort(vap.begin(), vap.end(), sort_tpd);
  std::sort(liq.begin(), liq.end(), sort_tpd);

  std::vector<Cand> selected;
  if (!vap.empty()) selected.push_back(vap.front());
  if (!liq.empty()) selected.push_back(liq.front());
  if (liq.size() >= 2) selected.push_back(liq[1]);

  // Fill remaining by global best TPD, keeping uniqueness in (phase_flag,x)
  for (const auto& u : uniq) {
    if (static_cast<int>(selected.size()) >= opt.max_incipient) break;
    bool exists = false;
    for (const auto& s : selected) {
      if (s.phase_flag == u.phase_flag && l1_dist(s.x, u.x) < opt.distinct_l1) {
        exists = true;
        break;
      }
    }
    if (!exists) selected.push_back(u);
  }

  if (static_cast<int>(selected.size()) > opt.max_incipient) {
    selected.resize(static_cast<size_t>(opt.max_incipient));
  }

  out.incipient.reserve(selected.size());
  for (const auto& s : selected) {
    IncipientPhase ip;
    ip.phase_flag = s.phase_flag;
    ip.reference_phase_flag = s.ref_flag;
    ip.x = s.x;
    ip.tpd = s.tpd;
    ip.iters = s.iters;
    out.incipient.push_back(std::move(ip));
  }

  if (opt.verbose) {
    std::cout << "[Stability] g_v=" << g_v << " g_l=" << g_l
              << " report_ref=" << out.reference_phase_flag
              << " unstable_candidates=" << uniq.size()
              << " returned=" << out.incipient.size() << "\n";
  }

  return out;

}

} // namespace phase_stability
