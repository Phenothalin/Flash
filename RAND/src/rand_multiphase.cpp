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
  // 1. 构建 Context
  SystemContext sys;
  sys.temperature = input.temperature;
  sys.pressure = input.pressure;
  sys.feedMoles = input.feedMoles;
  sys.elementMatrix = elementMatrix;
  sys.phases.resize(nPhases);

  // 2. 初始化逻辑
  if (nPhases == 2 && initialPhaseCompositions.empty()) {
      // 两相且无初值：使用 Wilson K-factor 估算 (调用 initializeTwoPhase)
      auto init = initializeTwoPhase(sys, {}, {});
      sys.phases[0].state = {input.temperature, input.pressure, init.n_phases[0], thermo_.vaporPhaseFlag()};
      sys.phases[1].state = {input.temperature, input.pressure, init.n_phases[1], thermo_.liquidPhaseFlag()};
  } 
  else {
      // === 【新增逻辑】多相初值估算：基于投影法估算 Beta ===
      // 目标：找到一组 beta，使得 sum(beta_j * x_j) ≈ feed
      // 投影近似公式：beta_j ≈ (x_j · feed) / (x_j · x_j)
      // 这种方法对于组成差异大的相（如水、油、气）非常有效且稳健
      
      std::vector<double> estimated_betas(nPhases, 0.0);
      double sum_beta_est = 0.0;
      double total_feed_moles = std::accumulate(input.feedMoles.begin(), input.feedMoles.end(), 0.0);

      for (int j = 0; j < nPhases; ++j) {
          if (j < (int)initialPhaseCompositions.size() && !initialPhaseCompositions[j].empty()) {
              const auto& x = initialPhaseCompositions[j];
              // 确保 x 是归一化的，方便计算
              double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
              std::vector<double> x_norm = x;
              if (sum_x > 1e-12) {
                  for (auto& val : x_norm) val /= sum_x;
              }

              double dot_feed_x = 0.0;
              double dot_x_x = 0.0;
              for (size_t i = 0; i < input.feedMoles.size(); ++i) {
                  dot_feed_x += input.feedMoles[i] * x_norm[i];
                  dot_x_x += x_norm[i] * x_norm[i];
              }

              // 投影估算 beta
              double beta = (dot_x_x > 1e-12) ? (dot_feed_x / dot_x_x) : 0.0;
              
              // 保护机制：防止 beta 过小或为负
              if (beta < 1e-4 * total_feed_moles) beta = 1e-4 * total_feed_moles;
              
              estimated_betas[j] = beta;
          } else {
              // 没有给初值的相，先给个极小值占位，或者均分
              estimated_betas[j] = total_feed_moles / nPhases; 
          }
          sum_beta_est += estimated_betas[j];
      }

      // 归一化 Betas 以严格匹配总进料量
      if (sum_beta_est > 1e-12) {
          double scale = total_feed_moles / sum_beta_est;
          for (auto& b : estimated_betas) b *= scale;
      }

      // 打印一下估算的 Beta，方便调试
      std::cout << "[Init] Estimated Betas: ";
      for(auto b : estimated_betas) std::cout << b << " ";
      std::cout << std::endl;

      // === 应用估算的 Beta 生成 n ===
      for(int j=0; j<nPhases; ++j) {
          std::vector<double> n;
          if (j < (int)initialPhaseCompositions.size() && !initialPhaseCompositions[j].empty()) {
              n = initialPhaseCompositions[j];
              // 归一化并乘以估算的 beta
              double s = std::accumulate(n.begin(), n.end(), 0.0);
              if (s > 1e-12) {
                  double factor = estimated_betas[j] / s;
                  for(auto& val : n) val *= factor;
              }
          } else {
              // Fallback: 均分进料 (理论上不应走到这里，如果初值都给了的话)
              n = input.feedMoles;
              double split = 1.0 / nPhases;
              for(auto& val : n) val *= split;
          }
          
          // 简单的相态标记分配
          int flag = (j==0 ? thermo_.vaporPhaseFlag() : thermo_.liquidPhaseFlag());
          sys.phases[j].state = {input.temperature, input.pressure, n, flag};
      }
  }

  // 3. 核心求解
  return solveGeneral(sys, maxIter, tol);
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

MultiFlashResult RandFlash::solveGeneral(
    SystemContext& sys,
    int maxIter,
    double tol)
{
    MultiFlashResult result;
    result.pressure = sys.pressure;
    result.temperature = sys.temperature;
    result.feedComposition = sys.feedMoles;
    
    const size_t F = sys.phases.size();
    const size_t E = sys.elementMatrix.size();

    std::vector<double> Acoef, rhs, sol;
    std::vector<double> Lambda;
    std::vector<double> deltaBeta(F);
    std::vector<std::vector<double>> dnPhases(F);

    std::cout << "进入通用迭代循环 (F=" << F << "), maxIter = " << maxIter << std::endl;

    try {
        for (int iter = 0; iter < maxIter; ++iter) {
            std::cout << "\n=== RandFlash General Iter " << (iter + 1) << " ===\n";

            std::vector<std::vector<std::vector<double>>> Ms(F);
            std::vector<std::vector<double>> mus(F);
            std::vector<std::vector<double>> nPhases(F);

            for(size_t j=0; j<F; ++j) {
                updatePhaseChemistry(sys.phases[j]);
                
                // 【核心修复】提高特征值下限阈值
                // 1e-10 -> 1e-3
                // 这限制了不稳定相的修正步长幅度，防止产生过大的 delta_n 导致 alpha 过小
                PhaseFixOptions pfx;
                pfx.eig_floor = 1e-3; 
                pfx.inv_tol   = 1e-12;
                
                fixPhaseHessian(sys.phases[j], pfx);

                Ms[j] = sys.phases[j].M;
                mus[j] = sys.phases[j].mu;
                nPhases[j] = sys.phases[j].state.moleNumbers;
            }

            assembleGlobalSystem(
                sys.temperature,
                Ms, mus, nPhases,
                sys.elementMatrix,
                sys.feedMoles,
                Acoef, rhs
            );

            double lin_resid = 0.0;
            sol = solveGlobalLinearSystem(Acoef, rhs, &lin_resid);
            
            // 抑制虚假的高残差警告（相对误差可能很小）
            // 只有当残差真的非常大时才打印
            // (feed - current) 产生的 RHS 项可能很大，导致 lin_resid 绝对值大，这是正常的
            if (lin_resid > 1e-4 * (1.0 + std::abs(sol[0]))) { // 简单的相对检查
                 // std::cerr << "[RandFlash] Linear resid: " << lin_resid << "\n";
            }

            if (sol.size() != E + F) throw std::runtime_error("Linear solve size mismatch");

            Lambda.assign(sol.begin(), sol.begin() + E);
            for(size_t j=0; j<F; ++j) deltaBeta[j] = sol[E + j];

            backSubstituteDeltas(
                sys.temperature,
                sys.elementMatrix,
                Ms, mus, nPhases,
                Lambda, deltaBeta,
                dnPhases
            );

            std::vector<std::vector<double>> currentNs(F);
            for(size_t j=0; j<F; ++j) currentNs[j] = sys.phases[j].state.moleNumbers;

            applyUpdate(sys.temperature, mus, dnPhases, currentNs);

            for(size_t j=0; j<F; ++j) sys.phases[j].state.moleNumbers = currentNs[j];

            auto conv = checkConvergence(
                iter + 1, tol, sys.temperature,
                mus, sys.elementMatrix, currentNs, sys.feedMoles
            );

            if (conv.converged) {
                result.success = true;
                result.iterations = iter + 1;
                result.mu_infinity_norm = conv.max_mu_diff;
                result.elem_residual_inf = conv.elem_error;
                
                result.n_phase = currentNs;
                result.beta.resize(F);
                for(size_t j=0; j<F; ++j) 
                    result.beta[j] = std::accumulate(currentNs[j].begin(), currentNs[j].end(), 0.0);
                
                return result;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "[RandFlash] Exception: " << e.what() << std::endl;
        result.success = false;
        return result;
    }

    result.success = false;
    return result;
}