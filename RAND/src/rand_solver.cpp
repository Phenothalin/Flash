// rand_plan.cpp
// 负责：相稳定性分析调用、相数/相型规划、求解入口 solve()（驱动初始化与通用求解器）、结果相序重排。
// 说明：该文件只做“流程编排”，数值内核在 rand_solver.cpp / rand_flash.cpp，初始化在 rand_init.cpp。

#include "rand_flash.hpp"
#include "phase_stability.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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
                        std::cout << "[PhaseMerge] merged phase " << kill
                                  << " into " << keep << ", L1(x)=" << dist << "\n";
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
            std::cout << "[PhasePrune] removed tiny phase " << kill
                      << " (beta=" << bmin << ") -> merged into " << tgt
                      << ", new F=" << sys.phases.size() << "\n";
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
    std::cout << "[Init] Initialized " << F << " phases with Least-Squares + Strict Element Scaling.\n";

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
        std::cout << "[Init] Initialized 1 phase with Reactive Element Conservation.\n";
    } else {
        init = initializeReactiveMultiPhase(sys, numPhases, phaseFlags, elementMoles);
        std::cout << "[Init] Initialized " << numPhases << " phases with Reactive Element Conservation.\n";
    }

    // 4) 填充相状态
    for (int j = 0; j < numPhases; ++j) {
        sys.phases[j].state = {input.temperature, input.pressure, init.n_phases[j], phaseFlags[j]};
    }

    // 5) 调用反应体系求解器
    auto res = solveGeneralReactive(sys, maxIter, tol);

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
// 通用接口：自动相稳定性分析 -> 自动相数选择 -> 自动初始化 -> 相分裂
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solve(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions,
    int maxIter,
    double tol)
{
    // 若用户给了初值，则按给定相数直接求解
    // 使用minGibbsPhaseFlag让ThermoPack自动选择Gibbs能最低的根
    if (!initialPhaseCompositions.empty()) {
        const int F = static_cast<int>(initialPhaseCompositions.size());
        std::vector<int> flags(F, thermo_.minGibbsPhaseFlag());
        return solveWithPhaseGuesses(input, elementMatrix, initialPhaseCompositions, flags, maxIter, tol);
    }

    // 1) 相稳定性分析
    phase_stability::PhaseStabilityAnalyzer analyzer(thermo_);
    phase_stability::StabilityOptions opt;
    opt.verbose = false;

    auto stab = analyzer.analyze(input.temperature, input.pressure, input.feedMoles, opt);

    // 2) 单相稳定：直接返回（不进入 Rand 迭代）
    if (stab.stable) {
        MultiFlashResult result;
        result.success = true;
        result.iterations = 0;
        result.pressure = input.pressure;
        result.temperature = input.temperature;
        result.feedComposition = input.feedMoles;
        result.mu_infinity_norm = 0.0;
        result.elem_residual_inf = 0.0;

        // 构造单相的PhaseContext，使用analyze中的参考相相态
        PhaseContext singlePhase;
        singlePhase.state.Temperature = input.temperature;
        singlePhase.state.Pressure = input.pressure;
        singlePhase.state.moleNumbers = input.feedMoles;
        singlePhase.state.phaseFlag = stab.reference_phase_flag;  // 使用稳定性分析的参考相
        result.phases = {singlePhase};

        return result;
    }

    // 3) 根据稳定性结果选择相数与相类型。
    // 关键约束：对“典型 VLE”体系（一个气相+一个液相）
    // 稳定性分析可能同时给出 vap 和 liq 的 incipient（数值噪声 / seed 多样性导致），
    // 但这不应被直接解释为 3 相。3 相仅在“需要两套液相”（LLE/LLV）时启用。

    auto normalize = [&](const std::vector<double>& v) {
        std::vector<double> x = v;
        for (double& xi : x) xi = std::max(xi, 1e-14);
        double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s <= 0.0) {
            const double uni = 1.0 / static_cast<double>(x.size());
            std::fill(x.begin(), x.end(), uni);
        } else {
            for (double& xi : x) xi /= s;
        }
        return x;
    };

    const int vap = thermo_.vaporPhaseFlag();
    const int liq = thermo_.liquidPhaseFlag();
    const int mingibbs = thermo_.minGibbsPhaseFlag();

    // Wilson K (用于在 stability 未给出某类 incipient 时构造更“有偏”的相组成初猜)
    std::vector<double> K(input.feedMoles.size(), 1.0);
    try {
        thermo_.wilsonK(input.temperature, input.pressure, K);
        for (double& Ki : K) Ki = std::clamp(Ki, 1e-12, 1e12);
    } catch (...) {
        // fallback: keep K=1
    }

    auto build_wilson_guess = [&](bool vapor_like) {
        std::vector<double> x = input.feedMoles;
        const double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s > 0.0) {
            for (double& xi : x) xi /= s;
        }
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = vapor_like ? (x[i] * K[i]) : (x[i] / K[i]);
        }
        return normalize(x);
    };

    // 探测是否含水（用于 LLV：水/油两液相）
    int water_idx = -1;
    {
        const auto names = thermo_.getComponentNames();
        for (size_t i = 0; i < std::min(names.size(), input.feedMoles.size()); ++i) {
            std::string n = names[i];
            std::transform(n.begin(), n.end(), n.begin(), ::toupper);
            if (n == "H2O" || n == "WATER" || n.find("H2O") != std::string::npos) {
                water_idx = static_cast<int>(i);
                break;
            }
        }
    }
    const double ztot = std::accumulate(input.feedMoles.begin(), input.feedMoles.end(), 0.0);
    const double z_water = (water_idx >= 0 && ztot > 0.0) ? (input.feedMoles[water_idx] / ztot) : 0.0;
    const bool has_water = (water_idx >= 0 && z_water > 1e-6);

    // 收集 incipient：最多 2 个
    std::vector<std::vector<double>> liquid_cands;
    std::vector<std::vector<double>> vapor_cands;
    for (const auto& inc : stab.incipient) {
        if (inc.phase_flag == liq) liquid_cands.push_back(inc.x);
        if (inc.phase_flag == vap) vapor_cands.push_back(inc.x);
    }

    // 默认：先做 2 相（最不容易引入回归）
    // VLE: {VAP, LIQ}
    // LLE: {LIQ, LIQ}
    // LLV: {VAP, LIQ, LIQ}

    // --- 2 相初猜 ---
    std::vector<std::vector<double>> comps2;
    std::vector<int> flags2;
    comps2.reserve(2);
    flags2.reserve(2);

    // vapor guess
    std::vector<double> xV = (!vapor_cands.empty()) ? vapor_cands.front() : build_wilson_guess(true);
    // liquid guess
    std::vector<double> xL = (!liquid_cands.empty()) ? liquid_cands.front() : build_wilson_guess(false);

    // 若 reference 是 vapor 或 liquid，则更倾向于把 z 放在对应相上（更稳健）
    if (stab.reference_phase_flag == vap) xV = normalize(input.feedMoles);
    if (stab.reference_phase_flag == liq) xL = normalize(input.feedMoles);

    // 判断是否更像 LLE（两套 liquid），仅当 stability 明确给出 >=2 个 liquid candidate
    const bool want_LLE = (liquid_cands.size() >= 2) && (vapor_cands.empty()) && !has_water;

    if (want_LLE) {
        comps2 = { liquid_cands[0], liquid_cands[1] };
        flags2 = { mingibbs, mingibbs };
    } else {
        comps2 = { xV, xL };
        flags2 = { mingibbs, mingibbs };
    }

    // --- 3 相判别（仅当确实需要两套 liquid 时才启用）---
    // 1) 水体系：优先尝试 LLV
    // 2) 非水体系：只有在 liquid candidate >= 2 时才尝试 LLV
    const bool want_LLV = has_water || (liquid_cands.size() >= 2);

    if (want_LLV && !want_LLE) {
        // 3 相初值：优先来自稳定性结果；不够则用通用 3 相初始化补齐
        SystemContext sys0;
        sys0.temperature = input.temperature;
        sys0.pressure = input.pressure;
        sys0.feedMoles = input.feedMoles;
        sys0.elementMatrix = elementMatrix;

        std::vector<std::vector<double>> comps3;
        std::vector<int> flags3 = {mingibbs, mingibbs, mingibbs};

        // vapor
        if (!vapor_cands.empty()) comps3.push_back(vapor_cands.front());
        else comps3.push_back(build_wilson_guess(true));

        // two liquids
        if (liquid_cands.size() >= 2) {
            comps3.push_back(liquid_cands[0]);
            comps3.push_back(liquid_cands[1]);
        } else {
            // 用 Generic 初始化构造两套液相种子（尤其对含水体系更稳健）
            auto init3 = initializeThreePhaseGeneric(sys0);
            // init3.compositions 是 {V, L1, L2}
            comps3 = init3.compositions;
        }

        auto res3 = solveWithPhaseGuesses(input, elementMatrix, comps3, flags3, maxIter, tol);
        if (res3.success) return res3;
        // 3 相失败时回退 2 相
    }

    return solveWithPhaseGuesses(input, elementMatrix, comps2, flags2, maxIter, tol);
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
    const int vap = thermo_.vaporPhaseFlag();
    const int liq = thermo_.liquidPhaseFlag();
    const int mingibbs = thermo_.minGibbsPhaseFlag();

    // 0) 若用户给了初始相组成：直接走给定相数（跳过 stability analysis）
    if (!opt.initial_phase_compositions.empty()) {
        const int F = static_cast<int>(opt.initial_phase_compositions.size());
        std::vector<int> flags;
        if (!opt.forced_phase_flags.empty()) {
            if (static_cast<int>(opt.forced_phase_flags.size()) != F) {
                throw std::invalid_argument("SolveOptions: forced_phase_flags size mismatch");
            }
            flags = opt.forced_phase_flags;
        } else {
            // 默认使用minGibbsPhaseFlag，让ThermoPack自动选择Gibbs能最低的根
            flags.assign(F, mingibbs);
        }
        return solveWithPhaseGuesses(input, elementMatrix, opt.initial_phase_compositions, flags, maxIter, tol);
    }

    // 1) enable_stability_test=true：保持原逻辑
    if (opt.enable_stability_test) {
        return this->solve(input, elementMatrix, std::vector<std::vector<double>>{}, maxIter, tol);
    }

    // 2) enable_stability_test=false：必须强制相数
    const int F = opt.forced_phase_count;
    if (F <= 0) {
        throw std::invalid_argument("SolveOptions: forced_phase_count must be > 0 when enable_stability_test=false and no initial compositions");
    }

    // 2.1 phase flags
    std::vector<int> flags;
    if (!opt.forced_phase_flags.empty()) {
        if (static_cast<int>(opt.forced_phase_flags.size()) != F) {
            throw std::invalid_argument("SolveOptions: forced_phase_flags size mismatch");
        }
        flags = opt.forced_phase_flags;
    } else {
        // 默认使用minGibbsPhaseFlag
        flags.assign(F, mingibbs);
    }

    // 2.2 helpers
    auto normalize = [&](const std::vector<double>& v) {
        std::vector<double> x = v;
        for (double& xi : x) xi = std::max(xi, 1e-14);
        double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s <= 0.0) {
            const double uni = 1.0 / static_cast<double>(x.size());
            std::fill(x.begin(), x.end(), uni);
        } else {
            for (double& xi : x) xi /= s;
        }
        return x;
    };

    // Wilson K (用于构造 vapor-like / liquid-like 初猜)
    std::vector<double> K(input.feedMoles.size(), 1.0);
    try {
        thermo_.wilsonK(input.temperature, input.pressure, K);
        for (double& Ki : K) Ki = std::clamp(Ki, 1e-12, 1e12);
    } catch (...) {
        // fallback: keep K=1
    }

    auto build_wilson_guess = [&](bool vapor_like) {
        std::vector<double> x = input.feedMoles;
        const double s = std::accumulate(x.begin(), x.end(), 0.0);
        if (s > 0.0) {
            for (double& xi : x) xi /= s;
        }
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = vapor_like ? (x[i] * K[i]) : (x[i] / K[i]);
        }
        return normalize(x);
    };

    auto perturb = [&](const std::vector<double>& base, int k) {
        std::vector<double> x = base;
        // deterministic small perturbation to break symmetry
        const double eps = 1e-3 * static_cast<double>(k + 1);
        for (size_t i = 0; i < x.size(); ++i) {
            const double sgn = (i % 2 == 0) ? 1.0 : -1.0;
            x[i] = std::max(1e-14, x[i] * (1.0 + eps * sgn));
        }
        return normalize(x);
    };

    const std::vector<double> z = normalize(input.feedMoles);
    const std::vector<double> xV = build_wilson_guess(true);
    const std::vector<double> xL = build_wilson_guess(false);

    // 2.3 generate phase composition guesses
    std::vector<std::vector<double>> comps;
    comps.reserve(F);

    if (F == 1) {
        // single-phase: use z as a neutral guess (phase flag decided by user / default)
        comps.push_back(z);
    } else if (F == 2) {
        for (int j = 0; j < F; ++j) {
            comps.push_back(flags[j] == vap ? xV : xL);
        }
    } else {
        // Use existing 3-phase generic initializer to create {V, L1, L2} seeds.
        SystemContext sys0;
        sys0.temperature = input.temperature;
        sys0.pressure = input.pressure;
        sys0.feedMoles = input.feedMoles;
        sys0.elementMatrix = elementMatrix;

        auto init3 = initializeThreePhaseGeneric(sys0);
        const std::vector<double> seedV  = (init3.compositions.size() >= 1) ? init3.compositions[0] : xV;
        const std::vector<double> seedL1 = (init3.compositions.size() >= 2) ? init3.compositions[1] : xL;
        const std::vector<double> seedL2 = (init3.compositions.size() >= 3) ? init3.compositions[2] : xL;

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

    return solveWithPhaseGuesses(input, elementMatrix, comps, flags, maxIter, tol);
}

MultiFlashResult RandFlash::solveGeneral(
    SystemContext& sys,
    int maxIter,
    double tol)
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

    std::cout << "进入通用迭代循环 (F=" << sys.phases.size() << "), maxIter = " << maxIter << std::endl;

    // Phase pruning/merge parameters (conservative defaults)
    const double beta_rel_prune = 1e-14;   // relative to total feed
    const double beta_abs_prune = 1e-16;   // absolute moles
    const double l1_merge_tol   = 1e-8;    // merge nearly identical phases

    try {
        for (int iter = 0; iter < maxIter; ++iter) {
            size_t F = sys.phases.size();
            if (F == 0) throw std::runtime_error("No phases in SystemContext");

            std::cout << "\n=== RandFlash General Iter " << (iter + 1) << " ===\n";

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

            for (size_t i = E; i < E + F; ++i) {
                std::cout << "∆β[" << (i - E) << "]=" << sol[i] << "  ";
            }
            std::cout << std::endl;

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

            // 5) Phase pruning/merge after update (handles beta->0 degeneracy)
            const bool changed = prune_and_merge_phases(sys, beta_rel_prune, beta_abs_prune, l1_merge_tol, /*verbose*/false);
            if (changed) {
                std::cout << "[PhasePrune] Active phases changed -> F=" << sys.phases.size() << std::endl;
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
                mus_post, sys.elementMatrix, n_post, sys.feedMoles);

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
    return result;
}

// ----------------------------------------------------------------------------------
// 反应体系求解内核：专用于反应体系的Newton迭代，使用步长收敛判据
// ----------------------------------------------------------------------------------
MultiFlashResult RandFlash::solveGeneralReactive(
    SystemContext& sys,
    int maxIter,
    double tol)
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

    std::cout << "进入反应体系迭代循环 (F=" << sys.phases.size() << "), maxIter = " << maxIter << std::endl;

    // Phase pruning/merge parameters (conservative defaults)
    const double beta_rel_prune = 1e-14;   // relative to total feed
    const double beta_abs_prune = 1e-16;   // absolute moles
    const double l1_merge_tol   = 1e-8;    // merge nearly identical phases

    try {
        for (int iter = 0; iter < maxIter; ++iter) {
            size_t F = sys.phases.size();
            if (F == 0) throw std::runtime_error("No phases in SystemContext");

            std::cout << "\n=== RandFlash Reactive Iter " << (iter + 1) << " ===\n";

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

            for (size_t i = E; i < E + F; ++i) {
                std::cout << "∆β[" << (i - E) << "]=" << sol[i] << "  ";
            }
            std::cout << std::endl;

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

            // 5) Phase pruning/merge after update (handles beta->0 degeneracy)
            const bool changed = prune_and_merge_phases(sys, beta_rel_prune, beta_abs_prune, l1_merge_tol, /*verbose*/false);
            if (changed) {
                std::cout << "[PhasePrune] Active phases changed -> F=" << sys.phases.size() << std::endl;
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

            // 【关键区别】传递 dnPhases 和 isReactive=true 给收敛判断
            auto conv = checkConvergence(
                iter + 1, tol, sys.temperature,
                mus_post, sys.elementMatrix, n_post, sys.feedMoles,
                dnPhases, true);  // Pass step size and reactive flag

            if (conv.converged) {
                result.success = true;
                result.iterations = iter + 1;
                result.mu_infinity_norm = conv.max_mu_diff;
                result.elem_residual_inf = conv.elem_error;

                result.phases = sys.phases;
                return result;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[RandFlash] Exception: " << e.what() << std::endl;
        result.success = false;
        result.phases = sys.phases;
        return result;
    }

    result.success = false;
    result.phases = sys.phases;
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

    std::cout << "\n=== Reactive System Auto Phase Number Determination ===\n";
    std::cout << "Max phases: " << maxPhases << "\n";
    std::cout << "Element moles: ";
    for (double em : elementMoles) std::cout << em << " ";
    std::cout << "\n\n";

    MultiFlashResult bestResult;
    double bestGibbs = 1e100;
    int bestPhaseCount = 1;

    const int mingibbs = thermo_.minGibbsPhaseFlag();

    // Try phase counts from 1 to maxPhases
    for (int F = 1; F <= maxPhases; ++F) {
        std::cout << "\n--- Trying F = " << F << " phases ---\n";

        // Setup solve options for forced phase count
        SolveOptions opt;
        opt.enable_stability_test = false;
        opt.forced_phase_count = F;
        opt.forced_phase_flags.assign(F, mingibbs);

        // Solve with F phases
        auto result = solve(input, elementMatrix, opt, maxIterations, tolerance);

        if (!result.success) {
            std::cout << "F = " << F << " failed to converge, stopping.\n";
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
        std::cout << "F = " << F << " converged, G_total = " << G << " J\n";

        // Check if this is better than previous
        if (F == 1 || G < bestGibbs - 1e-6) {
            // Significant improvement
            bestGibbs = G;
            bestResult = result;
            bestPhaseCount = F;
            std::cout << "  -> New best solution (ΔG = " << (F > 1 ? bestGibbs - G : 0.0) << ")\n";
        } else {
            // No improvement, stop
            std::cout << "  -> No improvement (ΔG = " << (G - bestGibbs) << "), stopping.\n";
            break;
        }
    }

    std::cout << "\n=== Auto Phase Determination Complete ===\n";
    std::cout << "Optimal phase count: " << bestPhaseCount << "\n";
    std::cout << "Final Gibbs energy: " << bestGibbs << " J\n\n";

    return bestResult;
}

} // namespace randflash
