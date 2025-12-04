#pragma once
#include <vector>
#include <memory>
#include <stdexcept>

namespace ls {

/// 通用线性求解接口：封装具体后端（Eigen / MKL / cuSolver 等）
/// 所有矩阵均按行主序展开：A[i*n + j] = 第 i 行第 j 列
struct LinearSolverInterface {
  virtual ~LinearSolverInterface() = default;

  /// 求解对称正定线性系统 A x = b
  virtual std::vector<double> solveSPD(
      int n,
      const std::vector<double>& A,
      const std::vector<double>& b) = 0;

  /// 求解一般稠密线性系统 A x = b（不要求对称 / 正定）
  /// 要求返回最小二乘意义下的“最佳”解；如果 residual_norm 非空，
  /// 则写入 ||A x - b||_2
  virtual std::vector<double> solveDense(
      int n,
      const std::vector<double>& A,
      const std::vector<double>& b,
      double* residual_norm = nullptr) = 0;

  /// 对称实矩阵特征分解
  ///  输入: A (n×n, 行主序)
  ///  输出: eigenvalues[0..n-1] 为升序特征值，
  ///        eigenvectors 为 n×n 行主序矩阵，列 j 为特征向量 v_j：
  ///        eigenvectors[i*n + j] = (v_j)_i
  virtual void eigenDecomposeSymmetric(
      int n,
      const std::vector<double>& A,
      std::vector<double>& eigenvalues,
      std::vector<double>& eigenvectors) = 0;

  /// 计算对称正定矩阵 A 的逆 A^{-1}
  ///  返回矩阵同样按行主序展开
  virtual std::vector<double> invertSPD(
      int n,
      const std::vector<double>& A) = 0;
};

/// 创建基于 Eigen 的默认线性求解器实现
std::unique_ptr<LinearSolverInterface> createEigenSolver();

/// ============ 一些简单的工具函数（仅使用 std::vector，不依赖任何后端） ============

/// 将二维 std::vector 矩阵按行主序压平
inline std::vector<double> flattenRowMajor(
    const std::vector<std::vector<double>>& M)
{
  const int n = static_cast<int>(M.size());
  if (n == 0) return {};
  const int m = static_cast<int>(M[0].size());
  std::vector<double> out;
  out.reserve(n * m);
  for (int i = 0; i < n; ++i) {
    if (static_cast<int>(M[i].size()) != m) {
      throw std::runtime_error("flattenRowMajor: non-rectangular matrix");
    }
    out.insert(out.end(), M[i].begin(), M[i].end());
  }
  return out;
}

/// 从行主序一维数组还原二维矩阵
inline std::vector<std::vector<double>> unflattenRowMajor(
    int n, int m, const std::vector<double>& data)
{
  if (static_cast<int>(data.size()) != n * m) {
    throw std::runtime_error("unflattenRowMajor: size mismatch");
  }
  std::vector<std::vector<double>> M(n, std::vector<double>(m));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      M[i][j] = data[i * m + j];
    }
  }
  return M;
}

/// 使用给定求解器对对称正定矩阵做求逆（二维矩阵版本）
/// 主要是给 RandFlash 这种使用 std::vector<std::vector<double>> 的代码调用
inline std::vector<std::vector<double>> invert(
    const std::vector<std::vector<double>>& M,
    LinearSolverInterface& solver)
{
  const int n = static_cast<int>(M.size());
  if (n == 0) return {};
  if (static_cast<int>(M[0].size()) != n) {
    throw std::runtime_error("invert: matrix must be square");
  }
  std::vector<double> A_flat = flattenRowMajor(M);
  std::vector<double> inv_flat = solver.invertSPD(n, A_flat);
  return unflattenRowMajor(n, n, inv_flat);
}

// === 构造切空间正交基 B: 1^T y = 0 ===
// 用 (e_i - e_C) 做初基，再 Gram–Schmidt 正交化
inline std::vector<std::vector<double>> tangentBasis(int C){
  std::vector<std::vector<double>> B(C, std::vector<double>(C-1, 0.0));

  // 初始基：第 k 列为 e_k - e_C
  for (int k = 0; k < C-1; ++k){ // 列 k
    B[k][k]   =  1.0;
    B[C-1][k] = -1.0;
  }

  // Gram-Schmidt 正交化
  for (int j = 0; j < C-1; ++j){
    // 去除在之前列上的分量
    for (int i = 0; i < j; ++i){
      double proj = 0.0;
      for (int r = 0; r < C; ++r) {
        proj += B[r][i] * B[r][j];
      }
      for (int r = 0; r < C; ++r) {
        B[r][j] -= proj * B[r][i];
      }
    }

    // 归一化
    double nrm2 = 0.0;
    for (int r = 0; r < C; ++r) {
      nrm2 += B[r][j] * B[r][j];
    }
    double nrm = std::sqrt(nrm2);
    if (nrm < 1e-14) { // 退化保护：随机扰动再正交
      for (int r = 0; r < C; ++r) {
        B[r][j] = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
      }
      for (int i = 0; i < j; ++i){
        double proj = 0.0;
        for (int r = 0; r < C; ++r) {
          proj += B[r][i] * B[r][j];
        }
        for (int r = 0; r < C; ++r) {
          B[r][j] -= proj * B[r][i];
        }
      }
      nrm2 = 0.0;
      for (int r = 0; r < C; ++r) {
        nrm2 += B[r][j] * B[r][j];
      }
      nrm = std::sqrt(nrm2);
      if (nrm < 1e-14) {
        throw std::runtime_error("tangentBasis: failed to build non-degenerate basis");
      }
    }
    const double inv_nrm = 1.0 / nrm;
    for (int r = 0; r < C; ++r) {
      B[r][j] *= inv_nrm;
    }
  }

  // 可选：检查 1^T B = 0，这里略
  return B;
}

} // namespace ls
