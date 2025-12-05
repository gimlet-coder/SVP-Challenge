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

void LLL(Matrix &B, const Scalar delta);

/**
 * @brief LLL基底簡約アルゴリズム
 * @param B [in, out] n次元格子 L の基底 B
 *  最後に delta に関する LLL簡約基底 B として返すため変更あり
 * @param delta [in] LLLにおける簡約パラメータ 変更無し
 */

void DeepLLL(Matrix &B, const Scalar delta);

/**
 * @brief DeepLLL基底簡約アルゴリズム
 * @param B [in, out] n次元格子 L の基底 B
 *  最後に delta に関する DeepLLL簡約基底 B として返すため変更あり
 * @param delta [in] DeepLLLにおける簡約パラメータ 変更無し
 */

void MLLL(Matrix &B, const Scalar delta);

/**
 * @brief MLLL基底簡約アルゴリズム
 * @param B [in, out] n次元格子 L の基底 B
 *  最後に delta に関する MLLL 簡約基底 B として返すため変更あり
 * @param delta [in] MLLLにおける簡約パラメータ 変更無し
 */

 bool ENUM(const Matrix &U, const Vector &B_norm, Scalar &R_square , Vector &v_out, const int k_begin, const int k_end, long long &node_count);

 /** 
  * @brief 格子上の最短ベクトルの数え上げ ENUM 
  *   返り値 
  *     true: R_square より短いベクトルが見つかった場合
  *     false 見つからなかった場合
  * @param U [in] GSO係数行列 mu_{i, j} が U(i - 1, j - 1) に格納されていることに注意 変更無し
  * @param B_norm [in] GSOベクトルの2乗ノルム 変更無し
  * @param R_square [in] 数え上げ上界列 R_n^2 = R_square と考える 関数内で更新していくので変更あり
  * @param v_out [out] 条件を満たす格子ベクトルの係数ベクトル 
  * @param k_begin, k_end [in] 探索範囲 変更無し
  * @param node_count [out] 最短ベクトルを探すのにかかったノード数を保存する
  */


void BKZ(Matrix &B, int beta, const Scalar delta);

/**
 * @brief BKZ基底簡約アルゴリズム
 * @param beta [in] ブロックサイズ (2 <= beta <= n) 変更無し
 * @param delta [in] LLLの簡約パラメータ (1/4 < delta < 1) 変更無し
 */

void Progressive_BKZ(const std::string& filename, int n, int start_beta, int max_beta, Scalar delta, double target_norm, int max_retries);

/**
 * @brief ProgressiveBKZアルゴリズム
 * @param filename [in] 読み込ませたい行列のファイル名を入れる 変更無し
 * @param n [in] 行列サイズ
 * @param start_beta, max_beta [in] それぞれ beta のスタート値と最大値を入力する
 * @param delta [in] LLL の簡約パラメータ
 * @param target_norm [in] 既に既知のデータを扱う場合、おおよその最小ノルムを入力してそこを目指して計算させる
 * @param max_retries [in] 繰り返し計算する上限回数
 */