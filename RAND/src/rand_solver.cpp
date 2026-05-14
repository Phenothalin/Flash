// rand_plan.cpp
// 负责：相稳定性分析调用、相数/相型规划、求解入口 solve()（驱动初始化与通用求解器）、结果相序重排。
// 说明：该文件只做"流程编排"，数值内核在 rand_solver.cpp / rand_flash.cpp，初始化在 rand_init.cpp。

#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include "rand_utils.hpp"
#include "phase_stability.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

namespace randflash {
namespace {
static int find_water_index(const thermo::IThermoBackend& thermo, size_t C) {
    const auto names = thermo.getComponentNames();
    for (size_t i = 0; i < std::min(names.size(), C); ++i) {
        std::string n = names[i];
        std::transform(n.begin(), n.end(), n.begin(), ::toupper);
        if (n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// 根据压缩因子确定实际相态并更新phaseFlag
// Z > 0.5 为气相，Z <= 0.5 为液相
static void determine_phase_types(
    MultiFlashResult& res,
    thermo::IThermoBackend& thermo)
{
    const int vap = thermo.vaporPhaseFlag();
    const int liq = thermo.liquidPhaseFlag();

    for (size_t j = 0; j < res.phases.size(); ++j) {
        double beta = res.beta(j);
        if (beta < 1e-10) continue;  // 跳过消失的相

        // 计算压缩因子
        double Z = thermo.compressibilityFactor(res.phases[j].state);

        // 根据Z值确定phaseFlag
        if (Z > 0.5) {
            res.phases[j].state.phaseFlag = vap;
        } else {
            res.phases[j].state.phaseFlag = liq;
        }
    }
}

static void reorder_flash_result(
    MultiFlashResult& res,
    int vaporFlag,
    int water_idx)
{
    const size_t F = res.phases.size();
    if (F <= 1) return;

    // Build (idx, isVap, beta, x_water)
    struct Key {
        size_t idx;
        bool isVap;
        double beta;
        double xw;
    };
    std::vector<Key> keys;
    keys.reserve(F);

    for (size_t j = 0; j < F; ++j) {
        const auto& n = res.phases[j].state.moleNumbers;
        double beta = std::accumulate(n.begin(), n.end(), 0.0);
        double xw = 0.0;
        if (water_idx >= 0 && beta > 0.0 && static_cast<size_t>(water_idx) < n.size()) {
            xw = n[static_cast<size_t>(water_idx)] / beta;
        }
        bool isVap = (res.phases[j].state.phaseFlag == vaporFlag);
        keys.push_back(Key{j, isVap, beta, xw});
    }

    // Primary sort: vapor first
    // Secondary: for liquids, if water exists and there are >=2 liquids, order by water fraction ascending
    // Otherwise: by beta descending
    const bool hasWater = (water_idx >= 0);
    const int nLiq = static_cast<int>(std::count_if(keys.begin(), keys.end(), [](const Key& k){ return !k.isVap; }));

    std::sort(keys.begin(), keys.end(), [&](const Key& a, const Key& b){
        if (a.isVap != b.isVap) return a.isVap > b.isVap;
        if (!a.isVap && !b.isVap && hasWater && nLiq >= 2) {
            // oil first (low xw), water-rich later
            if (std::abs(a.xw - b.xw) > 1e-8) return a.xw < b.xw;
        }
        return a.beta > b.beta;
    });

    // Apply permutation
    std::vector<PhaseContext> phases_new;
    phases_new.reserve(F);

    for (const auto& k : keys) {
        phases_new.push_back(res.phases[k.idx]);
    }
    res.phases.swap(phases_new);
}


static double phase_beta(const std::vector<double>& n) {
    return std::accumulate(n.begin(), n.end(), 0.0);
}

static double l1_dist_composition(const std::vector<double>& na, double beta_a,
                                 const std::vector<double>& nb, double beta_b) {
    if (beta_a <= 0.0 || beta_b <= 0.0) return 1.0e100;
    const size_t C = std::min(na.size(), nb.size());
    double s = 0.0;
    for (size_t i = 0; i < C; ++i) {
        const double xa = na[i] / beta_a;
        const double xb = nb[i] / beta_b;
        s += std::abs(xa - xb);
    }
    return s;
}

// 相消失/相合并（Phase pruning / merge）
// - 目的：当某相 beta→0 时，避免进入病态多相系统并产生“伪均匀相”。
// - 策略：
//   1) 先合并同相型且组成非常接近的相（避免重复相）；
//   2) 再裁剪 beta 很小的相，并把其摩尔数并入最合适的目标相。
static bool prune_and_merge_phases(
    randflash::SystemContext& sys,
    double beta_rel_prune,
    double beta_abs_prune,
    double l1_merge_tol,
    bool verbose)
{
    size_t F = sys.phases.size();
    if (F <= 1) return false;

    const size_t C = sys.feedMoles.size();
    double tot = std::accumulate(sys.feedMoles.begin(), sys.feedMoles.end(), 0.0);
    if (tot <= 0.0) tot = 1.0;

    bool changed = false;

    // ---- 1) Merge similar phases (same phaseFlag) ----
    for (;;) {
        F = sys.phases.size();
        bool merged = false;
        for (size_t i = 0; i < F && !merged; ++i) {
            for (size_t j = i + 1; j < F; ++j) {
                const int fi = sys.phases[i].state.phaseFlag;
                const int fj = sys.phases[j].state.phaseFlag;
                if (fi != fj) continue;

                const auto& ni = sys.phases[i].state.moleNumbers;
                const auto& nj = sys.phases[j].state.moleNumbers;
                const double bi = phase_beta(ni);
                const double bj = phase_beta(nj);

                const double dist = l1_dist_composition(ni, bi, nj, bj);
                if (dist < l1_merge_tol) {
                    const size_t keep = (bi >= bj) ? i : j;
                    const size_t kill = (bi >= bj) ? j : i;

                    for (size_t k = 0; k < C; ++k) {
                        sys.phases[keep].state.moleNumbers[k] += sys.phases[kill].state.moleNumbers[k];
                    }
                    sys.phases.erase(sys.phases.begin() + static_cast<long>(kill));
                    merged = true;
                    changed = true;

                    if (verbose) {
                        RAND_DEBUG("[PhaseMerge] merged phase {} into {}, L1(x)={}", kill, keep, dist);
                    }
                    break;
                }
            }
        }
        if (!merged) break;
    }

    // ---- 2) Prune tiny phases ----
    const double prune_threshold = std::max(beta_abs_prune, beta_rel_prune * tot);

    for (;;) {
        F = sys.phases.size();
        if (F <= 1) break;

        size_t kill = F;
        double bmin = 1.0e300;
        for (size_t j = 0; j < F; ++j) {
            const double bj = phase_beta(sys.phases[j].state.moleNumbers);
            if (bj < prune_threshold && bj < bmin) {
                bmin = bj;
                kill = j;
            }
        }
        if (kill == F) break;

        const int fkill = sys.phases[kill].state.phaseFlag;

        // Prefer merging into same phaseFlag with largest beta
        size_t tgt = F;
        double bbest = -1.0;
        for (size_t j = 0; j < F; ++j) {
            if (j == kill) continue;
            const double bj = phase_beta(sys.phases[j].state.moleNumbers);
            if (sys.phases[j].state.phaseFlag == fkill && bj > bbest) {
                bbest = bj;
                tgt = j;
            }
        }
        // Otherwise merge into overall largest beta
        if (tgt == F) {
            for (size_t j = 0; j < F; ++j) {
                if (j == kill) continue;
                const double bj = phase_beta(sys.phases[j].state.moleNumbers);
                if (bj > bbest) {
                    bbest = bj;
                    tgt = j;
                }
            }
        }

        if (tgt == F) break; // should not happen

        for (size_t k = 0; k < C; ++k) {
            sys.phases[tgt].state.moleNumbers[k] += sys.phases[kill].state.moleNumbers[k];
        }
        sys.phases.erase(sys.phases.begin() + static_cast<long>(kill));
        changed = true;

        if (verbose) {
            RAND_DEBUG("[PhasePrune] removed tiny phase {} (beta={}) -> merged into {}, new F={}",
                      kill, bmin, tgt, sys.phases.size());
        }
    }

    return changed;
}
} // namespace

// ----------------------------------------------------------------------------------
// 内部驱动：给定相数组合猜测 + 相标志，完成初始化 + solveGeneral + 结果相序稳定化
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solveWithPhaseGuesses(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& phaseCompositions,
    const std::vector<int>& phaseFlags,
    int maxIter,
    double tol)
{
    const int F = static_cast<int>(phaseCompositions.size());
    if (F <= 0) {
        throw std::invalid_argument("solveWithPhaseGuesses: phaseCompositions is empty");
    }
    if (static_cast<int>(phaseFlags.size()) != F) {
        throw std::invalid_argument("solveWithPhaseGuesses: phaseFlags size mismatch");
    }

    // 1) 构建系统上下文
    SystemContext sys;
    sys.temperature = input.temperature;
    sys.pressure = input.pressure;
    sys.feedMoles = input.feedMoles;
    sys.elementMatrix = elementMatrix;
    sys.phases.resize(F);

    // 2) 非反应体系初始化：使用物质守恒初始化
    InitResult init = initializeFromCompositions(sys, F, phaseCompositions, phaseFlags);
    RAND_INFO("[Init] Initialized {} phases with Least-Squares + Strict Element Scaling.", F);

    for (int j = 0; j < F; ++j) {
        sys.phases[j].state = {input.temperature, input.pressure, init.n_phases[j], phaseFlags[j]};
    }

    // 3) 核心求解
    auto res = solveGeneral(sys, maxIter, tol);

    // 4) 根据压缩因子确定实际相态
    if (res.success) {
        determine_phase_types(res, thermo_);
    }

    // 5) 结果相序稳定化（Vapor -> (Oil-like) -> (Water-like)）
    if (res.success && res.phases.size() == sys.phases.size()) {
        const int vap = thermo_.vaporPhaseFlag();
        const int water_idx = find_water_index(thermo_, input.feedMoles.size());
        reorder_flash_result(res, vap, water_idx);
    }
    return res;
}

// ----------------------------------------------------------------------------------
// 反应体系求解接口：独立的反应体系求解入口
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solveReactive(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    int numPhases,
    int maxIter,
    double tol)
{
    const size_t C = input.feedMoles.size();
    const size_t E = elementMatrix.size();

    if (numPhases < 1) {
        throw std::invalid_argument("solveReactive: numPhases must be at least 1");
    }

    // 1) 计算元素摩尔数
    std::vector<double> elementMoles(E, 0.0);
    for (size_t e = 0; e < E; ++e) {
        for (size_t i = 0; i < C; ++i) {
            elementMoles[e] += elementMatrix[e][i] * input.feedMoles[i];
        }
    }

    // 2) 构建系统上下文
    SystemContext sys;
    sys.temperature = input.temperature;
    sys.pressure = input.pressure;
    sys.feedMoles = input.feedMoles;
    sys.elementMatrix = elementMatrix;
    sys.phases.resize(numPhases);

    // 3) 反应体系初始化
    std::vector<int> phaseFlags(numPhases, thermo_.minGibbsPhaseFlag());
    InitResult init;

    if (numPhases == 1) {
        init = initializeReactiveSinglePhase(sys, elementMoles, phaseFlags[0]);
        RAND_INFO("[Init] Initialized 1 phase with Reactive Element Conservation.");
    } else {
        init = initializeReactiveMultiPhase(sys, numPhases, phaseFlags, elementMoles);
        RAND_INFO("[Init] Initialized {} phases with Reactive Element Conservation.", numPhases);
    }

    // 4) 填充相状态
    for (int j = 0; j < numPhases; ++j) {
        sys.phases[j].state = {input.temperature, input.pressure, init.n_phases[j], phaseFlags[j]};
    }

    // 5) 调用反应体系求解器
    auto res = solveGeneral(sys, maxIter, tol, true);  // isReactive=true

    // 6) 根据压缩因子确定实际相态
    if (res.success) {
        determine_phase_types(res, thermo_);
    }

    // 7) 结果相序稳定化（Vapor -> (Oil-like) -> (Water-like)）
    if (res.success && res.phases.size() == sys.phases.size()) {
        const int vap = thermo_.vaporPhaseFlag();
        const int water_idx = find_water_index(thermo_, input.feedMoles.size());
        reorder_flash_result(res, vap, water_idx);
    }

    return res;
}


// ----------------------------------------------------------------------------------
// 新接口：SolveOptions
// - 默认保持原 solve() 行为（enable_stability_test=true 且不提供初猜时：自动 stability analysis）
// - 支持关闭 stability analysis 并强制相数（reaction 体系常用）
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solve(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIter,
    double tol)
{
    // Pure routing logic: dispatch to reactive or non-reactive paths
    if (opt.is_reactive) {
        return solveReactiveDispatch(input, elementMatrix, opt, maxIter, tol);
    }
    return solveNonReactiveDispatch(input, elementMatrix, opt, maxIter, tol);
}

MultiFlashResult RandFlash::solveGeneral(
    SystemContext& sys,
    int maxIter,
    double tol,
    bool isReactive)
{
    MultiFlashResult result;
    result.pressure = sys.pressure;
    result.temperature = sys.temperature;
    result.feedComposition = sys.feedMoles;

    const size_t E = sys.elementMatrix.size();
    const size_t C = sys.feedMoles.size();

    std::vector<double> Acoef, rhs, sol;
    std::vector<double> Lambda;
    std::vector<double> deltaBeta;
    std::vector<std::vector<double>> dnPhases;

    RAND_DEBUG("进入通用迭代循环 (F={}), maxIter = {}", sys.phases.size(), maxIter);

    // Phase pruning/merge parameters (conservative defaults)
    const double beta_rel_prune = 1e-14;   // relative to total feed
    const double beta_abs_prune = 1e-16;   // absolute moles
    const double l1_merge_tol   = 1e-8;    // merge nearly identical phases

    try {
        for (int iter = 0; iter < maxIter; ++iter) {
            size_t F = sys.phases.size();
            if (F == 0) throw std::runtime_error("No phases in SystemContext");

            RAND_DEBUG("\n=== RandFlash General Iter {} ===", iter + 1);

            // Allocate/resize per-iteration containers
            deltaBeta.assign(F, 0.0);
            dnPhases.assign(F, std::vector<double>(C, 0.0));

            std::vector<std::vector<std::vector<double>>> Ms(F);
            std::vector<std::vector<double>> mus(F);
            std::vector<std::vector<double>> nPhases(F);

            // 1) Update chemistry & Hessian fix for each phase
            for (size_t j = 0; j < F; ++j) {
                updatePhaseChemistry(sys.phases[j]);

                PhaseFixOptions pfx;
                // 维持与两相版本一致的默认量级，并让 line-search 决定实际步长。
                pfx.eig_floor = 1e-10;
                pfx.inv_tol   = 1e-12;
                fixPhaseHessian(sys.phases[j], pfx);

                Ms[j]     = sys.phases[j].M;
                mus[j]    = sys.phases[j].mu;
                nPhases[j]= sys.phases[j].state.moleNumbers;
            }

            // 2) Assemble & solve global linear system
            assembleGlobalSystem(
                sys.temperature,
                Ms, mus, nPhases,
                sys.elementMatrix,
                sys.feedMoles,
                Acoef, rhs);

            double lin_resid = 0.0;
            sol = solveGlobalLinearSystem(Acoef, rhs, &lin_resid);

            std::string deltaBetaStr;
            for (size_t i = E; i < E + F; ++i) {
                deltaBetaStr += fmt::format("∆β[{}]={:.6g}  ", i - E, sol[i]);
            }
            RAND_DEBUG("{}", deltaBetaStr);

            if (sol.size() != E + F) throw std::runtime_error("Linear solve size mismatch");

            Lambda.assign(sol.begin(), sol.begin() + static_cast<long>(E));
            for (size_t j = 0; j < F; ++j) deltaBeta[j] = sol[E + j];

            // 3) Back-substitute to obtain dn for each phase
            backSubstituteDeltas(
                sys.temperature,
                sys.elementMatrix,
                Ms, mus, nPhases,
                Lambda, deltaBeta,
                dnPhases);

            // 4) Apply update with line search and constraints
            std::vector<std::vector<double>> currentNs(F);
            for (size_t j = 0; j < F; ++j) currentNs[j] = sys.phases[j].state.moleNumbers;

            applyUpdate(sys.temperature, mus, dnPhases, currentNs);
            for (size_t j = 0; j < F; ++j) sys.phases[j].state.moleNumbers = currentNs[j];

            // 4.5) 基于误差大小的元素守恒投影（统一策略）
            if (E > 0) {
                // Compute target element moles from feed
                std::vector<double> targetElementMoles(E, 0.0);
                for (size_t e = 0; e < E; ++e) {
                    for (size_t i = 0; i < C; ++i) {
                        targetElementMoles[e] += sys.elementMatrix[e][i] * sys.feedMoles[i];
                    }
                }

                // Compute current element conservation error
                std::vector<double> currentElementMoles(E, 0.0);
                for (size_t e = 0; e < E; ++e) {
                    for (size_t j = 0; j < F; ++j) {
                        for (size_t i = 0; i < C; ++i) {
                            currentElementMoles[e] += sys.elementMatrix[e][i] * currentNs[j][i];
                        }
                    }
                }

                double elem_err = 0.0;
                for (size_t e = 0; e < E; ++e) {
                    elem_err = std::max(elem_err, std::abs(currentElementMoles[e] - targetElementMoles[e]));
                }

                // Apply projection only if error exceeds threshold
                // Use 1e-6 as threshold to match convergence criterion and avoid
                // creating numerical noise that prevents convergence
                const double projection_threshold = 1e-6;
                if (elem_err > projection_threshold) {
                    RAND_DEBUG("  [Conservation] Element error {:.6e} > threshold, applying projection", elem_err);
                    projectToElementConservation(sys.elementMatrix, targetElementMoles, currentNs);
                    for (size_t j = 0; j < F; ++j) sys.phases[j].state.moleNumbers = currentNs[j];
                }
            }

            // 5) Phase pruning/merge after update (handles beta->0 degeneracy)
            // 对于反应体系，禁用相合并以保持元素守恒
            if (!isReactive) {
                const bool changed = prune_and_merge_phases(sys, beta_rel_prune, beta_abs_prune, l1_merge_tol, /*verbose*/false);
                if (changed) {
                    RAND_DEBUG("[PhasePrune] Active phases changed -> F={}", sys.phases.size());
                }
            }

            // 6) Recompute mu for convergence check using the UPDATED state
            F = sys.phases.size();
            std::vector<std::vector<double>> mus_post(F);
            std::vector<std::vector<double>> n_post(F);
            for (size_t j = 0; j < F; ++j) {
                updatePhaseChemistry(sys.phases[j]);
                mus_post[j] = sys.phases[j].mu;
                n_post[j]   = sys.phases[j].state.moleNumbers;
            }

            auto conv = checkConvergence(
                iter + 1, tol, sys.temperature,
                mus_post, sys.elementMatrix, n_post, sys.feedMoles,
                dnPhases, isReactive, Lambda);

            result.iter_mu_history.push_back(conv.max_mu_diff);

            if (conv.converged) {
                result.success = true;
                result.iterations = iter + 1;
                result.mu_infinity_norm = conv.max_mu_diff;
                result.elem_residual_inf = conv.elem_error;

                result.phases = sys.phases;  // 保存各相的完整上下文（含phaseFlag）
                return result;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[RandFlash] Exception: " << e.what() << std::endl;
        result.success = false;
        result.phases = sys.phases;  // 保存当前状态
        return result;
    }

    result.success = false;
    result.phases = sys.phases;  // 保存当前状态

    // Check for local minimum trap in reactive single-phase systems
    if (isReactive && sys.phases.size() == 1 && !sys.elementMatrix.empty()) {
        if (isLikelyLocalMinimum(sys, sys.elementMatrix)) {
            RAND_INFO("WARNING: Solver appears trapped at local minimum (uniform distribution).");
            RAND_INFO("This occurs because PR/SRK EOS lacks standard Gibbs formation energies.");
            RAND_INFO("Consider using a thermodynamic model with ΔG_f° for reactive equilibrium.");
        }
    }

    return result;
}

// ----------------------------------------------------------------------------------
// 反应体系自动相数判断接口：顺序相添加法
// 从单相开始，逐步尝试添加新相，直到Gibbs能不再降低
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solveReactiveAuto(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    int maxPhases,
    int maxIterations,
    double tolerance)
{
    if (maxPhases < 1) {
        throw std::invalid_argument("maxPhases must be at least 1");
    }

    const size_t C = input.feedMoles.size();
    const size_t E = elementMatrix.size();

    // Compute element moles from feed composition
    std::vector<double> elementMoles(E, 0.0);
    for (size_t e = 0; e < E; ++e) {
        for (size_t i = 0; i < C; ++i) {
            elementMoles[e] += elementMatrix[e][i] * input.feedMoles[i];
        }
    }

    std::string elemStr;
    for (double em : elementMoles) elemStr += fmt::format("{:.6g} ", em);
    RAND_INFO("\n=== Reactive System Auto Phase Number Determination ===");
    RAND_INFO("Max phases: {}", maxPhases);
    RAND_INFO("Element moles: {}", elemStr);

    MultiFlashResult bestResult;
    double bestGibbs = 1e100;
    int bestPhaseCount = 1;

    const int mingibbs = thermo_.minGibbsPhaseFlag();

    // Try phase counts from 1 to maxPhases
    for (int F = 1; F <= maxPhases; ++F) {
        RAND_DEBUG("\n--- Trying F = {} phases ---", F);

        // Setup solve options for forced phase count
        SolveOptions opt;
        opt.enable_stability_test = false;
        opt.forced_phase_count = F;
        opt.forced_phase_flags.assign(F, mingibbs);

        // Solve with F phases
        auto result = solve(input, elementMatrix, opt, maxIterations, tolerance);

        if (!result.success) {
            RAND_DEBUG("F = {} failed to converge, stopping.", F);
            break;
        }

        // Build SystemContext to compute Gibbs energy
        SystemContext sys;
        sys.temperature = input.temperature;
        sys.pressure = input.pressure;
        sys.feedMoles = input.feedMoles;
        sys.elementMatrix = elementMatrix;
        sys.phases = result.phases;

        // Ensure chemical potentials are computed
        for (size_t j = 0; j < sys.phases.size(); ++j) {
            updatePhaseChemistry(sys.phases[j]);
        }

        double G = computeTotalGibbs(sys);
        RAND_DEBUG("F = {} converged, G_total = {} J", F, G);

        // Check if this is better than previous
        if (F == 1 || G < bestGibbs - 1e-6) {
            // Significant improvement
            bestGibbs = G;
            bestResult = result;
            bestPhaseCount = F;
            RAND_DEBUG("  -> New best solution (ΔG = {})", (F > 1 ? bestGibbs - G : 0.0));
        } else {
            // No improvement, stop
            RAND_DEBUG("  -> No improvement (ΔG = {}), stopping.", (G - bestGibbs));
            break;
        }
    }

    RAND_INFO("\n=== Auto Phase Determination Complete ===");
    RAND_INFO("Optimal phase count: {}", bestPhaseCount);
    RAND_INFO("Final Gibbs energy: {} J", bestGibbs);

    return bestResult;
}

// ----------------------------------------------------------------------------------
// Architecture Refactoring: Routing Functions
// ----------------------------------------------------------------------------------

MultiFlashResult RandFlash::solveReactiveDispatch(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations,
    double tolerance)
{
    if (opt.auto_phase_count || opt.forced_phase_count <= 0) {
        return solveReactiveAuto(input, elementMatrix, opt.max_phases, maxIterations, tolerance);
    }
    return solveReactive(input, elementMatrix, opt.forced_phase_count, maxIterations, tolerance);
}

MultiFlashResult RandFlash::solveNonReactiveDispatch(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations,
    double tolerance)
{
    // Determine if we should use automatic phase determination
    bool auto_phase = opt.enable_stability_test &&
                      opt.initial_phase_compositions.empty() &&
                      opt.forced_phase_count <= 0;

    if (auto_phase) {
        return solveNonReactiveAuto(input, elementMatrix, opt, maxIterations, tolerance);
    }
    return solveNonReactiveFixed(input, elementMatrix, opt, maxIterations, tolerance);
}

MultiFlashResult RandFlash::solveNonReactiveAuto(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations,
    double tolerance)
{
    if (opt.phase_determination_strategy == "stable") {
        return solveStable(input, elementMatrix, opt, maxIterations, tolerance);
    }
    return solveFast(input, elementMatrix, opt, maxIterations, tolerance);
}

MultiFlashResult RandFlash::solveNonReactiveFixed(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations,
    double tolerance)
{
    // User provided initial compositions
    if (!opt.initial_phase_compositions.empty()) {
        const int F = static_cast<int>(opt.initial_phase_compositions.size());
        std::vector<int> flags;

        if (!opt.forced_phase_flags.empty()) {
            if (static_cast<int>(opt.forced_phase_flags.size()) != F) {
                throw std::invalid_argument("SolveOptions: forced_phase_flags size mismatch");
            }
            flags = opt.forced_phase_flags;
        } else {
            // Default: use minGibbsPhaseFlag
            flags.assign(F, thermo_.minGibbsPhaseFlag());
        }

        return solveWithPhaseGuesses(input, elementMatrix,
            opt.initial_phase_compositions, flags, maxIterations, tolerance);
    }

    // No initial compositions, use default initialization
    return solveWithDefaultInit(input, elementMatrix, opt, maxIterations, tolerance);
}

MultiFlashResult RandFlash::solveWithDefaultInit(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations,
    double tolerance)
{
    const int F = opt.forced_phase_count;
    if (F <= 0) {
        throw std::invalid_argument(
            "SolveOptions: forced_phase_count must be > 0 when enable_stability_test=false and no initial compositions");
    }

    const int vap = thermo_.vaporPhaseFlag();
    const int liq = thermo_.liquidPhaseFlag();
    const int mingibbs = thermo_.minGibbsPhaseFlag();

    // Determine phase flags
    std::vector<int> flags;
    if (!opt.forced_phase_flags.empty()) {
        if (static_cast<int>(opt.forced_phase_flags.size()) != F) {
            throw std::invalid_argument("SolveOptions: forced_phase_flags size mismatch");
        }
        flags = opt.forced_phase_flags;
    } else {
        // Default: use minGibbsPhaseFlag
        flags.assign(F, mingibbs);
    }

    // Get Wilson K-values for composition guesses
    std::vector<double> K(input.feedMoles.size(), 1.0);
    try {
        thermo_.wilsonK(input.temperature, input.pressure, K);
        for (double& Ki : K) Ki = std::clamp(Ki, 1e-12, 1e12);
    } catch (...) {
        // fallback: keep K=1
    }

    const std::vector<double> z = utils::normalize(input.feedMoles);
    const std::vector<double> xV = utils::buildWilsonGuess(input.feedMoles, K, true);
    const std::vector<double> xL = utils::buildWilsonGuess(input.feedMoles, K, false);

    // Generate phase composition guesses
    std::vector<std::vector<double>> comps;
    comps.reserve(F);

    if (F == 1) {
        // Single-phase: use z as neutral guess
        comps.push_back(z);
    } else if (F == 2) {
        // Two-phase: use vapor-like and liquid-like guesses
        for (int j = 0; j < F; ++j) {
            comps.push_back(flags[j] == vap ? xV : xL);
        }
    } else {
        // Multi-phase: use three-phase generic initializer
        SystemContext sys0;
        sys0.temperature = input.temperature;
        sys0.pressure = input.pressure;
        sys0.feedMoles = input.feedMoles;
        sys0.elementMatrix = elementMatrix;

        auto init3 = initializeThreePhaseGeneric(sys0);
        const std::vector<double> seedV  = (init3.compositions.size() >= 1) ? init3.compositions[0] : xV;
        const std::vector<double> seedL1 = (init3.compositions.size() >= 2) ? init3.compositions[1] : xL;
        const std::vector<double> seedL2 = (init3.compositions.size() >= 3) ? init3.compositions[2] : xL;

        // Perturb function for breaking symmetry
        auto perturb = [](const std::vector<double>& base, int k) {
            std::vector<double> x = base;
            const double eps = 1e-3 * static_cast<double>(k + 1);
            for (size_t i = 0; i < x.size(); ++i) {
                const double sgn = (i % 2 == 0) ? 1.0 : -1.0;
                x[i] = std::max(1e-14, x[i] * (1.0 + eps * sgn));
            }
            return utils::normalize(x);
        };

        int vap_used = 0;
        int liq_used = 0;
        const int liq_count = static_cast<int>(std::count(flags.begin(), flags.end(), liq));

        for (int j = 0; j < F; ++j) {
            if (flags[j] == vap) {
                if (vap_used == 0) comps.push_back(seedV);
                else comps.push_back(perturb(seedV, vap_used));
                ++vap_used;
            } else {
                if (liq_count >= 2) {
                    if (liq_used == 0) comps.push_back(seedL1);
                    else if (liq_used == 1) comps.push_back(seedL2);
                    else comps.push_back(perturb(seedL2, liq_used - 1));
                } else {
                    comps.push_back(xL);
                }
                ++liq_used;
            }
        }
    }

    return solveWithPhaseGuesses(input, elementMatrix, comps, flags, maxIterations, tolerance);
}

MultiFlashResult RandFlash::solveFast(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const SolveOptions& opt,
    int maxIterations,
    double tolerance)
{
    // Fast strategy: single stability analysis + direct phase splitting
    // This is the original solve() v1 logic

    // 1) Phase stability analysis
    phase_stability::PhaseStabilityAnalyzer analyzer(thermo_);
    phase_stability::StabilityOptions stabOpt;
    stabOpt.verbose = false;

    auto stab = analyzer.analyze(input.temperature, input.pressure, input.feedMoles, stabOpt);
    for(const auto& inc : stab.incipient) {
        RAND_DEBUG("[Stability] Found incipient phase: flag={}, TPD={}, iters={}",
                  inc.phase_flag, inc.tpd, inc.iters);
    }

    // 2) Single-phase stable: return immediately
    if (stab.stable) {
        MultiFlashResult result;
        result.success = true;
        result.iterations = 0;
        result.pressure = input.pressure;
        result.temperature = input.temperature;
        result.feedComposition = input.feedMoles;
        result.mu_infinity_norm = 0.0;
        result.elem_residual_inf = 0.0;

        PhaseContext singlePhase;
        singlePhase.state.Temperature = input.temperature;
        singlePhase.state.Pressure = input.pressure;
        singlePhase.state.moleNumbers = input.feedMoles;
        singlePhase.state.phaseFlag = stab.reference_phase_flag;
        result.phases = {singlePhase};

        return result;
    }

    // 3) Multi-phase: scenario classification based on stability results
    const int vap = thermo_.vaporPhaseFlag();
    const int liq = thermo_.liquidPhaseFlag();

    // Get Wilson K-values for fallback guesses
    std::vector<double> K(input.feedMoles.size(), 1.0);
    try {
        thermo_.wilsonK(input.temperature, input.pressure, K);
        for (double& Ki : K) Ki = std::clamp(Ki, 1e-12, 1e12);
    } catch (...) {
        // fallback: keep K=1
    }

    // Collect incipient phases
    std::vector<std::vector<double>> liquid_cands;
    std::vector<std::vector<double>> vapor_cands;
    for (const auto& inc : stab.incipient) {
        if (inc.phase_flag == liq) liquid_cands.push_back(inc.x);
        if (inc.phase_flag == vap) vapor_cands.push_back(inc.x);
    }

    // Vapor and liquid guesses (with fallback to Wilson)
    std::vector<double> xV = (!vapor_cands.empty()) ?
        vapor_cands.front() : utils::buildWilsonGuess(input.feedMoles, K, true);
    std::vector<double> xL = (!liquid_cands.empty()) ?
        liquid_cands.front() : utils::buildWilsonGuess(input.feedMoles, K, false);

    // If reference phase is vapor or liquid, use feed as that phase
    if (stab.reference_phase_flag == vap) xV = utils::normalize(input.feedMoles);
    if (stab.reference_phase_flag == liq) xL = utils::normalize(input.feedMoles);

    // Scenario classification
    const bool is_VLLE = (!vapor_cands.empty()) && (liquid_cands.size() >= 2);
    const bool is_LLE = (vapor_cands.empty()) && (liquid_cands.size() >= 2);

    // Scenario 1: VLLE - try 3-phase first
    if (is_VLLE) {
        std::vector<std::vector<double>> comps3;
        std::vector<int> flags3 = {vap, liq, liq};

        comps3.push_back(vapor_cands.front());
        comps3.push_back(liquid_cands[0]);
        comps3.push_back(liquid_cands[1]);

        auto res3 = solveWithPhaseGuesses(input, elementMatrix, comps3, flags3, maxIterations, tolerance);
        if (res3.success) return res3;

        // 3-phase failed, fallback to 2-phase VLE
        std::vector<std::vector<double>> comps2 = { xV, xL };
        std::vector<int> flags2 = { vap, liq };
        return solveWithPhaseGuesses(input, elementMatrix, comps2, flags2, maxIterations, tolerance);
    }

    // Scenario 2: LLE - 2-phase liquid-liquid
    if (is_LLE) {
        std::vector<std::vector<double>> comps2 = { liquid_cands[0], liquid_cands[1] };
        std::vector<int> flags2 = { liq, liq };
        return solveWithPhaseGuesses(input, elementMatrix, comps2, flags2, maxIterations, tolerance);
    }

    // Scenario 3: VLE - default 2-phase vapor-liquid
    std::vector<std::vector<double>> comps2 = { xV, xL };
    std::vector<int> flags2 = { vap, liq };
    return solveWithPhaseGuesses(input, elementMatrix, comps2, flags2, maxIterations, tolerance);
}

} // namespace randflash
