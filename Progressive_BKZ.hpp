#pragma once

// --- 必要なヘッダーファイル ---
// Eigenのコア機能とDenseモジュール（VectorXd, MatrixXdなど）を使うために必要
#include <Eigen/Dense>
#include "lattice_types.hpp"

// --- 関数の宣言 (プロトタイプ) ---



void DeepBKZ(Matrix &B, int beta, const Scalar delta);

/**
 * @brief DeepBKZ基底簡約アルゴリズム
 * @param beta [in] ブロックサイズ (2 <= beta <= n) 変更無し
 * @param delta [in] LLLの簡約パラメータ (1/4 < delta < 1) 変更無し
 */