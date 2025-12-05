#pragma once

// --- 必要なヘッダーファイル ---
// Eigenのコア機能とDenseモジュール（VectorXd, MatrixXdなど）を使うために必要
#include <Eigen/Dense>
#include "lattice_types.hpp"
// --- 関数の宣言 (プロトタイプ) ---



void Gram_Schmidt(const Matrix& B, Matrix& B_star, Matrix& U);
/**
 * @brief Gram-Schmidt の直交化を実行し、直交基底 B_star と係数行列 U を計算する
 * @param B [in] 入力となる基底行列 (n x m) 各行が基底ベクトル 変更無し
 * @param B_star [out] 計算された直交基底ベクトルを格納する行列 (n x m)。各行が b_i^*
 * @param U [out] GSO係数 μ_{i,j} (i > j) と対角成分 1 を格納する下三角行列 (n x n),
 *  μ_{i, j} が U(i - 1, j - 1) (i > j) に格納されることに注意 (0-based index)
 */

void LLL(Matrix &B, const Scalar delta);

/**
 * @brief LLL基底簡約アルゴリズム
 * @param B [in, out] n次元格子 L の基底 B
 *  最後に delta に関する LLL簡約基底 B として返すため変更あり
 * @param delta [in] LLLにおける簡約パラメータ 変更無し
 */

 void MLLL(Matrix &B, const Scalar delta);

/**
 * @brief MLLL基底簡約アルゴリズム
 * @param B [in, out] n次元格子 L の基底 B
 *  最後に delta に関する MLLL 簡約基底 B として返すため変更あり
 * @param delta [in] MLLLにおける簡約パラメータ 変更無し
 */

 bool ENUM(const Matrix &U, const Vector &B_norm, const Vector &R_squares , Vector &v_out, const int k_begin, const int k_end, long long &node_count);

 /** 
  * @brief 格子上の最短ベクトルの数え上げ ENUM 
  *   返り値 
  *     true: R_square より短いベクトルが見つかった場合
  *     false 見つからなかった場合
  * @param U [in] GSO係数行列 mu_{i, j} が U(i - 1, j - 1) に格納されていることに注意 変更無し
  * @param B_norm [in] GSOベクトルの2乗ノルム 変更無し
  * @param R_squares [in] 数え上げ上界列 枝刈りをして小さくしていく 変更無し
  * @param v_out [out] 条件を満たす格子ベクトルの係数ベクトル 
  * @param k_begin, k_end [in] 探索範囲 変更無し
  * @param node_count [out] 最短ベクトルを探すのにかかったノード数を保存する
  */