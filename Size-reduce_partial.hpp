#pragma once

// --- 必要なヘッダーファイル ---
// Eigenのコア機能とDenseモジュール（VectorXd, MatrixXdなど）を使うために必要
#include <Eigen/Dense>
#include "lattice_types.hpp"

// --- 関数の宣言 (プロトタイプ) ---

void Gram_Schmidt(const IntMatrix& B, RealMatrix& B_star, RealMatrix& U);
/**
 * @brief Gram-Schmidt の直交化を実行し、直交基底 B_star と係数行列 U を計算する
 * @param B [in] 入力となる基底行列 (n x m) 各行が基底ベクトル 変更無し
 * @param B_star [out] 計算された直交基底ベクトルを格納する行列 (n x m)。各行が b_i^*
 * @param U [out] GSO係数 μ_{i,j} (i > j) と対角成分 1 を格納する下三角行列 (n x n),
 *  μ_{i, j} が U(i - 1, j - 1) (i > j) に格納されることに注意 (0-based index)
 */


