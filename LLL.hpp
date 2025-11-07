#pragma once

// --- 必要なヘッダーファイル ---
// Eigenのコア機能とDenseモジュール（VectorXd, MatrixXdなど）を使うために必要
#include <Eigen/Dense>

using Vector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;

// --- 関数の宣言 (プロトタイプ) ---

void Gram_Schmidt(const Matrix& B, Matrix& B_star, Matrix& U);
/**
 * @brief Gram-Schmidt の直交化を実行し、直交基底 B_star と係数行列 U を計算する
 * @param B [in] 入力となる基底行列 (n x m) 各行が基底ベクトル 変更無し
 * @param B_star [out] 計算された直交基底ベクトルを格納する行列 (n x m)。各行が b_i^*
 * @param U [out] GSO係数 μ_{i,j} (i > j) と対角成分 1 を格納する下三角行列 (n x n),
 *  μ_{i, j} が U(i - 1, j - 1) (i > j) に格納されることに注意 (0-based index)
 */

void Size_reduce_partial(Matrix &B, Matrix &U, const int i,const int j);

/**
 * @brief アルゴリズム3に基づいて (i, j) 要素 において部分的にサイズ基底簡約をする
 * @param B [in, out] 入力となる基底行列 (n x m) 各行が基底ベクトル 関数内で変更する可能性あり
 * @param U [in, out] GSO係数 μ_{k, l} (k > l) と対角成分 1 を格納する下三角行列 (n x n),
 *  μ{k, l} が U(k - 1, l - 1) (k > l) が格納されることに注意 (0-based index) 関数内で変更する可能性あり
 * @param i, j [in] 更新対象の行インデックス (0-based) i > j である必要あり 変更無し
 */

void GSOUpdate_LLL_partial(Matrix &U, Vector &B_norm, const int k);

 /**
  * @brief アルゴリズム6に基づいて k における LLL内のGSO情報を更新する
  * @param U [in, out] GSO係数行列 mu_{i, j} が U(i - 1, j - 1) に格納されていることに注意
  *  更新されたGSO係数行列を U に反映させることから関数内で変更する可能性あり
  * @param B_norm [in, out] GSOベクトルの2乗ノルム B_i = ||b_i*||^2 が B_norm(i - 1) に格納されていることに注意 (0-based index)
  *  更新されたGSOベクトルの2乗ノルムを B_norm に反映させることから関数内で変更される可能性あり
  * @param k [in] 交換させるインデックス 2 <= k <= n に注意  変更無し
  */