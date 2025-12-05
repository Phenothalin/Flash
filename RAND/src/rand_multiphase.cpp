#include "rand_flash.hpp"
#include <numeric>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <stdexcept>
using namespace randflash;
using namespace ls;

MultiFlashResult RandFlash::solveMultiPhaseCore(
    const FlashInput& input,
    int nPhases,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions,
    int maxIter,
    double tol)
{
  MultiFlashResult out;
  out.pressure        = input.pressure;
  out.temperature     = input.temperature;
  out.feedComposition = input.feedMoles;

  // 目前核心算法只实现了两相，所以这里先硬性限制 nPhases == 2
  if (nPhases != 2) {
    std::cerr << "[RandFlash] solveMultiPhaseCore currently supports only 2 phases, "
              << "but nPhases = " << nPhases << std::endl;
    out.success = false;
    return out;
  }

  // 拆出两相初始组成（如果有的话），否则传空给 solveTwoPhase
  std::vector<double> y0, x0;
  if (!initialPhaseCompositions.empty()) {
    if (initialPhaseCompositions.size() >= 1)
      y0 = initialPhaseCompositions[0];
    if (initialPhaseCompositions.size() >= 2)
      x0 = initialPhaseCompositions[1];
  }

  // 调用现有两相求解器（Newton 核心仍在 solveTwoPhase 里）
  FlashResult fr = solveTwoPhase(
      input,
      elementMatrix,
      y0,
      x0,
      maxIter,
      tol);

  // 直接把失败状态返回
  if (!fr.success) {
    out.success           = false;
    out.iterations        = static_cast<int>(fr.iterations);
    out.mu_infinity_norm  = fr.convergenceError; // 这里是 max_mu_diff(dimless)
    out.elem_residual_inf = 0.0;                 // 下面只在成功时计算
    return out;
  }

  // === 成功：将两相结果填入 MultiFlashResult ===
  out.success          = true;
  out.iterations       = static_cast<int>(fr.iterations);
  out.mu_infinity_norm = fr.convergenceError;

  // n_phase: 绝对摩尔数（solveTwoPhase 在收敛时 result.vaporComposition / liquidComposition
  // 中存的就是 nV / nL，而不是归一化 x）
  out.n_phase.clear();
  out.n_phase.resize(2);
  out.n_phase[0] = fr.vaporComposition;   // n^V
  out.n_phase[1] = fr.liquidComposition;  // n^L

  // beta: 每相总摩尔数
  out.beta.resize(2);
  out.beta[0] = std::accumulate(out.n_phase[0].begin(), out.n_phase[0].end(), 0.0);
  out.beta[1] = std::accumulate(out.n_phase[1].begin(), out.n_phase[1].end(), 0.0);

  // 元素残差范数 ||A*(sum_j n^j - feed)||_inf
  if (!elementMatrix.empty()) {
    const size_t E = elementMatrix.size();
    const size_t C = elementMatrix[0].size();
    std::vector<double> n_sum(C, 0.0);

    // sum_j n^j
    for (int j = 0; j < 2; ++j) {
      for (size_t i = 0; i < C; ++i) {
        n_sum[i] += out.n_phase[j][i];
      }
    }

    // r_elem = A * (sum n - feed)
    double max_abs = 0.0;
    for (size_t ell = 0; ell < E; ++ell) {
      double val = 0.0;
      for (size_t i = 0; i < C; ++i) {
        val += elementMatrix[ell][i]
             * (n_sum[i] - input.feedMoles[i]);
      }
      max_abs = std::max(max_abs, std::fabs(val));
    }
    out.elem_residual_inf = max_abs;
  } else {
    out.elem_residual_inf = 0.0;
  }

  return out;
}

MultiFlashResult RandFlash::solveMultiPhase(
    const FlashInput& input,
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<std::vector<double>>& initialPhaseCompositions,
    int maxIter,
    double tol)
{
  // 目前：如果用户没给初值，就默认 2 相；
  // 如果给了，就用给的相数；但核心仍只接受 2 相。
  int F = initialPhaseCompositions.empty()
        ? 2
        : static_cast<int>(initialPhaseCompositions.size());

  return solveMultiPhaseCore(
      input,
      F,
      elementMatrix,
      initialPhaseCompositions,
      maxIter,
      tol);
}

