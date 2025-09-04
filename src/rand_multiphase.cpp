
// ============================== rand_multiphase.cpp ===============================
// New implementation file. Add it to your build alongside rand_flash.cpp.
// This file depends on rand_flash.hpp (for declarations, R_CONST, enums) and
// thermo_adapter.hpp (for PropertyPackageAdapter).

#include "rand_flash.hpp"
#include "thermo_adapter.hpp"
#include <Eigen/Dense>
#include <numeric>
#include <optional>
#include <cassert>
#include <iostream>
#include <limits>
#include <cmath>

namespace randflash {

// --------------------------- RandMultiphase::Impl ----------------------------
struct RandMultiphase::Impl {
    // Constructor wiring
    Impl(PropertyPackageType pkgType,
         std::shared_ptr<material_object::Cluster> cluster,
         ls::LinearSolverInterface& lin,
         const Options& opt)
    : pkgType_(pkgType), cluster_(std::move(cluster)), linearSolver_(lin), opt_(opt) {}

    // Public entry used by wrapper
    MultiFlashResult solvePT(
        double pressure,
        double temperature,
        const std::vector<double>& feed,                            // size C (ABS moles)
        const std::vector<std::vector<double>>& A,                  // E x C
        const std::vector<std::vector<double>>& initial_x)          // optional
    {
        const size_t C = feed.size();
        const double Ftot = sum_(feed);
        if (Ftot <= 0.0) throw std::invalid_argument("Feed total must be positive");

        // 1) init phases
        std::vector<std::vector<double>> n;   // F x C
        std::vector<double> beta;             // F
        initPhases_(feed, initial_x, n, beta);
        ensurePhaseModels_(n.size());

        if (opt_.verbose) {
            std::cout << "[RandMultiphase] init F=" << n.size() << ", beta fractions: ";
            double sumB = sum_(beta);
            for (double b : beta) std::cout << (b/sumB) << ' ';
            std::cout << '\n';
        }

        // 2) outer cycles: allow add/remove
        for (int cycle = 0; cycle < opt_.max_addremove_cycles; ++cycle) {
            auto [ok, n_solved, beta_solved, mu_spread, elem_inf] =
                newtonFixedF_(pressure, temperature, feed, A, n, beta);

            n = std::move(n_solved);
            beta = std::move(beta_solved);

            // 3) remove tiny phases
            const double remove_abs = opt_.remove_phase_tol_rel * Ftot;
            bool removed = false;
            for (size_t j = 0; j < n.size(); ) {
                if (beta[j] < remove_abs) {
                    if (opt_.verbose) std::cout << "  Remove phase j="<<j<<" beta="<<beta[j]<<"\n";
                    n.erase(n.begin()+j);
                    beta.erase(beta.begin()+j);
                    removed = true;
                } else {
                    ++j;
                }
            }
            if (removed) { ensurePhaseModels_(n.size()); continue; }

            // 4) stability hook: try to add 1 phase (minimal TPD)
            if (n.size() < static_cast<size_t>(opt_.max_phases)) {
                auto cand = proposeNewPhaseByStability_(pressure, temperature, n);
                if (cand.has_value()) {
                    std::vector<double> x = cand.value();
                    normalize_(x);
                    const double epsB = std::max(1e-6 * Ftot, 1e-9);
                    std::vector<double> n_new(x.size());
                    for (size_t i = 0; i < x.size(); ++i) n_new[i] = epsB * x[i];
                    n.push_back(std::move(n_new));
                    beta.push_back(epsB);
                    ensurePhaseModels_(n.size());
                    if (opt_.verbose) std::cout << "  Added phase from TPD (beta="<<epsB<<")\n";
                    continue; // re-solve with F+1
                }
            }

            // 5) done
            MultiFlashResult res;
            res.success = ok;
            res.pressure = pressure;
            res.temperature = temperature;
            res.feedComposition = feed;
            res.n_phase = std::move(n);
            res.beta = std::move(beta);
            res.iterations = last_newton_iters_;
            res.mu_infinity_norm = mu_spread;
            res.elem_residual_inf = elem_inf;
            return res;
        }

        MultiFlashResult res;
        res.success = false;
        res.pressure = pressure;
        res.temperature = temperature;
        res.feedComposition = feed;
        res.n_phase = std::move(n);
        res.beta = std::move(beta);
        res.iterations = last_newton_iters_;
        return res;
    }

    // ------------------------ core Newton (fixed F) ------------------------
    std::tuple<bool, std::vector<std::vector<double>>, std::vector<double>, double, double>
    newtonFixedF_(
        double P, double T,
        const std::vector<double>& feed,
        const std::vector<std::vector<double>>& A,
        std::vector<std::vector<double>> n,
        std::vector<double> beta)
    {
        const size_t C = feed.size();
        const size_t E = A.size();
        size_t F = n.size();

        std::vector<std::vector<double>> m_j, M_j; // CxC
        std::vector<std::vector<std::vector<double>>> M(F); // F x C x C
        std::vector<std::vector<double>> mu(F);

        std::vector<double> Acoef; // (E+F)^2 row-major
        std::vector<double> rhs;   // (E+F)
        std::vector<double> sol;   // (E+F)
        std::vector<double> Lambda;   // E
        std::vector<double> dBeta;    // F
        std::vector<std::vector<double>> dn; // F x C

        const double RT = R_CONST * T;
        last_newton_iters_ = 0;

        auto elemResidualInf = [&](const std::vector<std::vector<double>>& A,
                                   const std::vector<std::vector<double>>& n,
                                   const std::vector<double>& feed) {
            const size_t E = A.size();
            const size_t C = feed.size();
            double infn = 0.0;
            for (size_t ell = 0; ell < E; ++ell) {
                double acc = 0.0;
                for (size_t i = 0; i < C; ++i) {
                    double nsum = 0.0; for (size_t j = 0; j < n.size(); ++j) nsum += n[j][i];
                    acc += A[ell][i] * (nsum - feed[i]);
                }
                infn = std::max(infn, std::abs(acc));
            }
            return infn;
        };

        for (int it = 0; it < opt_.max_newton_iter; ++it) {
            // 1) local Jacobians & μ
            for (size_t j = 0; j < F; ++j) {
                assembleLocalJacobian_(T, P, n[j], phaseModels_[j], m_j, mu[j]);
                M[j] = invert_(m_j);
                beta[j] = sum_(n[j]);
            }

            // 2) global system (E+F)x(E+F)
            assembleGlobalSystemMulti_(T, M, mu, n, A, Acoef, rhs);

            const size_t N = E + F;
            assert(Acoef.size() == N*N && rhs.size() == N);
            Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
                A_mat(Acoef.data(), N, N);
            Eigen::Map<const Eigen::VectorXd> b_vec(rhs.data(), N);

            Eigen::MatrixXd A_reg = A_mat;
            double traceA = A_reg.trace();
            double jitter = std::max(1e-12, opt_.jitter_scale * std::abs(traceA));
            A_reg.diagonal().array() += jitter;

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_reg, Eigen::ComputeThinU | Eigen::ComputeThinV);
            Eigen::VectorXd x_vec = svd.solve(b_vec);

            sol.assign(x_vec.data(), x_vec.data() + N);
            Lambda.assign(sol.begin(), sol.begin() + E);
            dBeta.assign(sol.begin() + E, sol.end());

            // 3) back-substitute Δn
            backSubstituteDeltas_(T, A, Lambda, dBeta, n, mu, M, dn);

            // 4) line search
            double alpha = lineSearchMulti_(n, dn, mu, /*init*/1.0, opt_.line_alpha_min, opt_.line_shrink);

            // 5) update
            for (size_t j = 0; j < F; ++j)
                for (size_t i = 0; i < C; ++i)
                    n[j][i] += alpha * dn[j][i];
            for (size_t j = 0; j < F; ++j) beta[j] = sum_(n[j]);

            // 6) convergence checks
            double mu_spread = 0.0;
            for (size_t i = 0; i < C; ++i) {
                double min_mu =  std::numeric_limits<double>::infinity();
                double max_mu = -std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < F; ++j) {
                    min_mu = std::min(min_mu, mu[j][i]);
                    max_mu = std::max(max_mu, mu[j][i]);
                }
                mu_spread = std::max(mu_spread, std::abs(max_mu - min_mu));
            }
            double elem_inf = elemResidualInf(A, n, feed);

            if (opt_.verbose) {
                std::cout << "  it=" << (it+1)
                          << " | mu_spread=" << mu_spread
                          << " | elem_inf=" << elem_inf
                          << " | alpha=" << alpha << "\n";
            }

            last_newton_iters_ = it + 1;
            if (mu_spread < opt_.tol_mu && elem_inf < opt_.tol_elem) {
                return {true, n, beta, mu_spread, elem_inf};
            }
        }

        // not converged
        double mu_spread = 0.0;
        for (size_t i = 0; i < C; ++i) {
            double min_mu =  std::numeric_limits<double>::infinity();
            double max_mu = -std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < F; ++j) {
                min_mu = std::min(min_mu, mu[j][i]);
                max_mu = std::max(max_mu, mu[j][i]);
            }
            mu_spread = std::max(mu_spread, std::abs(max_mu - min_mu));
        }
        double elem_inf = 0.0; // compute with current n
        {
            const size_t E2 = A.size();
            const size_t C2 = feed.size();
            double infn = 0.0;
            for (size_t ell = 0; ell < E2; ++ell) {
                double acc = 0.0;
                for (size_t i = 0; i < C2; ++i) {
                    double nsum = 0.0; for (size_t j = 0; j < n.size(); ++j) nsum += n[j][i];
                    acc += A[ell][i] * (nsum - feed[i]);
                }
                infn = std::max(infn, std::abs(acc));
            }
            elem_inf = infn;
        }
        return {false, n, beta, mu_spread, elem_inf};
    }

    // --------------------------- Stability (TPD) ----------------------------
    // Minimal Michelsen-style TPD using only μ from the adapter
    // For each existing phase j with composition x_ref, define
    //   tau_i = μ_i(x_ref)/(RT) - ln x_ref_i
    // Iterate for trial composition w using fixed-point:
    //   w_i^{new} ∝ exp( tau_i - μ_i(w)/RT )
    // If min TPD < -tpd_eps -> return the corresponding w as candidate new phase.
    std::optional<std::vector<double>> proposeNewPhaseByStability_(
        double P, double T,
        const std::vector<std::vector<double>>& n_current)
    {
        const size_t F = n_current.size();
        if (F == 0) return std::nullopt;
        const size_t C = n_current[0].size();
        const double RT = R_CONST * T;

        double best_tpd = 0.0;
        std::optional<std::vector<double>> best_w;

        // Model for evaluating μ at the trial (new) phase
        thermo::PropertyPackageAdapter newModel(pkgType_, cluster_);
        newModel.setPhaseIndex(static_cast<int>(F)); // fresh index

        for (size_t j = 0; j < F; ++j) {
            // x_ref from current phase j
            std::vector<double> x_ref = comp_(n_current[j]);
            floorPositive_(x_ref, 1e-14); normalize_(x_ref);

            // μ at reference
            std::vector<double> mu_ref;
            {
                std::vector<double> n_ref = n_current[j];
                floorPositive_(n_ref, 1e-14);
                thermo::PhaseState st{ T, P, n_ref };
                mu_ref = phaseModels_[j].chemicalPotentials(st);
            }

            // tau_i = μ_ref/RT - ln x_ref
            std::vector<double> tau(C, 0.0);
            for (size_t i = 0; i < C; ++i) tau[i] = mu_ref[i]/RT - std::log(x_ref[i]);

            // initial w: one FP step from x_ref
            std::vector<double> w = x_ref;
            {
                std::vector<double> mu_w;
                std::vector<double> n_trial(C); for (size_t i = 0; i < C; ++i) n_trial[i] = w[i];
                thermo::PhaseState st{ T, P, n_trial };
                mu_w = newModel.chemicalPotentials(st);
                for (size_t i = 0; i < C; ++i) w[i] = std::exp(tau[i] - mu_w[i]/RT);
                normalize_(w);
            }

            // FP iterations
            const int    max_it   = 60;
            const double tol_w    = 1e-12;
            const double tpd_eps  = 1e-9;  // instability threshold
            for (int it = 0; it < max_it; ++it) {
                std::vector<double> n_trial(C); for (size_t i = 0; i < C; ++i) n_trial[i] = w[i];
                thermo::PhaseState st{ T, P, n_trial };
                std::vector<double> mu_w = newModel.chemicalPotentials(st);

                std::vector<double> w_new(C);
                for (size_t i = 0; i < C; ++i) w_new[i] = std::exp(tau[i] - mu_w[i]/RT);
                normalize_(w_new);

                double diff = maxAbsDiff_(w_new, w);
                w.swap(w_new);
                if (diff < tol_w) break;
            }

            // compute TPD(w) = sum_i w_i [ (μ_i(w)/RT - ln w_i) - tau_i ]
            double tpd = 0.0;
            {
                std::vector<double> n_trial(C); for (size_t i = 0; i < C; ++i) n_trial[i] = w[i];
                thermo::PhaseState st{ T, P, n_trial };
                std::vector<double> mu_w = newModel.chemicalPotentials(st);
                for (size_t i = 0; i < C; ++i) {
                    double gamma_w = mu_w[i]/RT - std::log(std::max(w[i], 1e-300));
                    tpd += w[i] * (gamma_w - tau[i]);
                }
            }

            if (j == 0 || tpd < best_tpd) {
                best_tpd = tpd;
                best_w = w;
            }
        }

        if (best_w.has_value() && best_tpd < -1e-9) return best_w;
        return std::nullopt;
    }

    // ------------------------------ Local ops ------------------------------
    void assembleLocalJacobian_(
        double T, double P,
        std::vector<double>& n_j,
        thermo::PropertyPackageAdapter& model,
        std::vector<std::vector<double>>& m_j,
        std::vector<double>& mu_j) const
    {
        thermo::PhaseState state{ T, P, n_j };
        mu_j = model.chemicalPotentials(state);
        floorPositive_(n_j, 1e-12);
        auto dmun = model.dMu_dN(state);
        const double beta = sum_(n_j);
        const double RT = R_CONST * T;
        const size_t C = n_j.size();
        m_j.assign(C, std::vector<double>(C, 0.0));
        for (size_t i = 0; i < C; ++i)
            for (size_t k = 0; k < C; ++k)
                m_j[i][k] = beta * (dmun[i][k] / RT) + 1.0;  // consistent with two-phase
    }

    void assembleGlobalSystemMulti_(
        double T,
        const std::vector<std::vector<std::vector<double>>>& M,    // F x C x C
        const std::vector<std::vector<double>>& mu,                // F x C
        const std::vector<std::vector<double>>& n,                 // F x C
        const std::vector<std::vector<double>>& A,                 // E x C
        std::vector<double>& Ac,                                   // row-major (E+F) x (E+F)
        std::vector<double>& rhs) const
    {
        const size_t F = n.size();
        const size_t E = A.size();
        const size_t C = n[0].size();
        const size_t N = E + F;
        const double RT = R_CONST * T;

        Ac.assign(N*N, 0.0);
        rhs.assign(N, 0.0);

        std::vector<double> beta(F, 0.0);
        std::vector<std::vector<double>> x(F, std::vector<double>(C, 0.0));
        for (size_t j = 0; j < F; ++j) {
            beta[j] = sum_(n[j]);
            for (size_t i = 0; i < C; ++i) x[j][i] = n[j][i] / beta[j];
        }

        // UL block: sum_j A (beta_j M_j) A^T
        for (size_t ell = 0; ell < E; ++ell) {
            for (size_t k = 0; k < E; ++k) {
                double s = 0.0;
                for (size_t j = 0; j < F; ++j) {
                    const auto& Mj = M[j];
                    const double b = beta[j];
                    for (size_t i = 0; i < C; ++i) {
                        double A_ell_i = A[ell][i]; if (A_ell_i == 0.0) continue;
                        for (size_t p = 0; p < C; ++p) {
                            double A_k_p = A[k][p]; if (A_k_p == 0.0) continue;
                            s += A_ell_i * (b * Mj[i][p]) * A_k_p;
                        }
                    }
                }
                Ac[ell*N + k] = s;
            }
        }

        // UR/LL blocks: A x_j
        for (size_t ell = 0; ell < E; ++ell) {
            for (size_t j = 0; j < F; ++j) {
                double s = 0.0;
                for (size_t i = 0; i < C; ++i) s += A[ell][i] * x[j][i];
                Ac[ell*N + (E + j)] = s;
                Ac[(E + j)*N + ell] = s;
            }
        }
        // LR block: zeros

        // RHS top: sum_j A (beta_j M_j (mu_j/RT))
        for (size_t ell = 0; ell < E; ++ell) {
            double val = 0.0;
            for (size_t j = 0; j < F; ++j) {
                const auto& Mj = M[j];
                const double b = beta[j];
                for (size_t i = 0; i < C; ++i) {
                    double A_ell_i = A[ell][i]; if (A_ell_i == 0.0) continue;
                    for (size_t p = 0; p < C; ++p)
                        val += A_ell_i * (b * Mj[i][p] * (mu[j][p] / RT));
                }
            }
            rhs[ell] = val;
        }

        // RHS bottom: red_j = x_j · (mu_j/RT)
        for (size_t j = 0; j < F; ++j) {
            double red = 0.0;
            for (size_t i = 0; i < C; ++i) red += x[j][i] * (mu[j][i] / RT);
            rhs[E + j] = red;
        }
    }

    void backSubstituteDeltas_(
        double T,
        const std::vector<std::vector<double>>& A,          // E x C
        const std::vector<double>& Lambda,                  // size E
        const std::vector<double>& deltaBeta,               // size F
        const std::vector<std::vector<double>>& n,          // F x C
        const std::vector<std::vector<double>>& mu,         // F x C
        const std::vector<std::vector<std::vector<double>>>& M, // F x C x C
        std::vector<std::vector<double>>& dn) const         // F x C
    {
        const size_t F = n.size();
        const size_t C = n[0].size();
        const size_t E = A.size();
        const double RT = R_CONST * T;

        std::vector<double> lambdaComp(C, 0.0);
        for (size_t p = 0; p < C; ++p)
            for (size_t ell = 0; ell < E; ++ell)
                lambdaComp[p] += A[ell][p] * Lambda[ell];

        dn.assign(F, std::vector<double>(C, 0.0));
        for (size_t j = 0; j < F; ++j) {
            const double beta = sum_(n[j]);
            std::vector<double> xj = comp_(n[j]);
            std::vector<double> diff(C, 0.0);
            for (size_t p = 0; p < C; ++p) diff[p] = lambdaComp[p] - (mu[j][p] / RT);
            for (size_t i = 0; i < C; ++i) {
                double comb = 0.0;
                for (size_t p = 0; p < C; ++p) comb += M[j][i][p] * diff[p];
                dn[j][i] = xj[i] * deltaBeta[j] + beta * comb;
            }
        }
    }

    double lineSearchMulti_(
        const std::vector<std::vector<double>>& n,
        const std::vector<std::vector<double>>& dn,
        const std::vector<std::vector<double>>& mu, // F x C
        double init_alpha, double min_alpha, double shrink) const
    {
        auto dot = [](const std::vector<double>& a, const std::vector<double>& b){
            double s = 0.0; for (size_t i = 0; i < a.size(); ++i) s += a[i]*b[i]; return s; };
        const size_t F = n.size();

        double dir = 0.0; // sum_j dn_j · μ_j  (RT factor only affects scale)
        for (size_t j = 0; j < F; ++j) dir += dot(dn[j], mu[j]);
        const double dir_eps = 1e-12;

        auto positive_after = [&](double a){
            for (size_t j = 0; j < F; ++j)
                for (size_t i = 0; i < n[j].size(); ++i)
                    if (n[j][i] + a*dn[j][i] <= 0.0) return false;
            return true;
        };

        double alpha = init_alpha;
        while (alpha > min_alpha) {
            bool positive = positive_after(alpha);
            bool descent  = (dir < -dir_eps) ? (alpha * dir < 0.0)
                                             : (std::abs(dir) <= dir_eps);
            if (positive && descent) break;
            alpha *= shrink;
        }
        if (alpha <= min_alpha) alpha = min_alpha;
        return alpha;
    }

    // ------------------------------- Helpers -------------------------------
    static void normalize_(std::vector<double>& x) {
        double s = sum_(x);
        if (!(s > 0.0)) throw std::invalid_argument("normalize_: sum <= 0");
        for (double& v : x) v /= s;
    }
    static double sum_(const std::vector<double>& v) {
        return std::accumulate(v.begin(), v.end(), 0.0);
    }
    static std::vector<double> comp_(const std::vector<double>& n) {
        const double b = sum_(n);
        std::vector<double> x(n.size(), 0.0);
        if (b > 0.0) for (size_t i = 0; i < n.size(); ++i) x[i] = n[i] / b;
        return x;
    }
    static void floorPositive_(std::vector<double>& v, double eps) {
        for (double& x : v) if (x < eps) x = eps;
    }
    static double maxAbsDiff_(const std::vector<double>& a, const std::vector<double>& b) {
        double m = 0.0; for (size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i]-b[i])); return m;
    }

    static std::vector<std::vector<double>> invert_(const std::vector<std::vector<double>>& mat) {
        const size_t n = mat.size();
        Eigen::MatrixXd M(n, n);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                M(i,j) = mat[i][j];
        double diag = M.diagonal().cwiseAbs().mean();
        double jitter = std::max(1e-12, 1e-10 * diag);
        M.diagonal().array() += jitter;
        Eigen::MatrixXd Mi = M.inverse();
        std::vector<std::vector<double>> out(n, std::vector<double>(n));
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                out[i][j] = Mi(i,j);
        return out;
    }

    void initPhases_(
        const std::vector<double>& feed,
        const std::vector<std::vector<double>>& initial_x,
        std::vector<std::vector<double>>& n_out,
        std::vector<double>& beta_out) const
    {
        const size_t C = feed.size();
        const double Ftot = sum_(feed);
        if (!initial_x.empty()) {
            const size_t F = initial_x.size();
            n_out.assign(F, std::vector<double>(C, 0.0));
            beta_out.assign(F, Ftot / static_cast<double>(F));
            for (size_t j = 0; j < F; ++j) {
                std::vector<double> x = initial_x[j];
                floorPositive_(x, 1e-14); normalize_(x);
                for (size_t i = 0; i < C; ++i) n_out[j][i] = beta_out[j] * x[i];
            }
            return;
        }
        // default: 2 phases, 50/50, x=z
        std::vector<double> z = feed; normalize_(z);
        const double beta = 0.5 * Ftot;
        n_out.assign(2, std::vector<double>(C, 0.0));
        beta_out = {beta, beta};
        for (size_t i = 0; i < C; ++i) {
            n_out[0][i] = beta * z[i];
            n_out[1][i] = beta * z[i];
        }
    }

    void ensurePhaseModels_(size_t F) {
        if (phaseModels_.size() == F) {
            for (size_t j = 0; j < F; ++j) phaseModels_[j].setPhaseIndex(static_cast<int>(j));
            return;
        }
        phaseModels_.clear();
        phaseModels_.reserve(F);
        for (size_t j = 0; j < F; ++j) {
            thermo::PropertyPackageAdapter model(pkgType_, cluster_);
            model.setPhaseIndex(static_cast<int>(j));
            phaseModels_.push_back(std::move(model));
        }
    }

    // Members
    PropertyPackageType pkgType_;
    std::shared_ptr<material_object::Cluster> cluster_;
    ls::LinearSolverInterface& linearSolver_;
    Options opt_;

    std::vector<thermo::PropertyPackageAdapter> phaseModels_;
    int last_newton_iters_ = 0;
};

// --------------------------- RandMultiphase wrapper ---------------------------
RandMultiphase::RandMultiphase(
    PropertyPackageType packageType,
    std::shared_ptr<material_object::Cluster> componentCluster,
    ls::LinearSolverInterface& linearSolver,
    const Options& opt)
: impl_(std::make_unique<Impl>(packageType, std::move(componentCluster), linearSolver, opt)) {}

RandMultiphase::~RandMultiphase() = default;

MultiFlashResult RandMultiphase::solvePT(
    double pressure,
    double temperature,
    const std::vector<double>& feed,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initial_x)
{
    return impl_->solvePT(pressure, temperature, feed, elementMatrix, initial_x);
}

} // namespace randflash
