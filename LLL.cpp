#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている



#include "LLL.hpp"
// LLL基底簡約で隣接するベクトルの交換後のGSO情報を更新する関数 (外部定義)
void GSOUpdate_LLL_partial(Matrix &U,Vector &B_norm, const int k); 
// Size-reduce_partial 関数 (外部定義)
void Size_reduce_partial(Matrix &B, Matrix &U, const int i_in, const int j_in);
// Gram_Schmidt 関数 (外部定義)
void Gram_Schmidt(const Matrix &B, Matrix &B_star, Matrix &U);


// アルゴリズム 5: LLL 基底簡約アルゴリズム
void LLL(Matrix &B, const Scalar delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size-reduce 不正なインデックス delta です (1/4 < delta < 1)");
    }
    // --- 引数チェック終了 ---

    int n = B.rows();
    if(n < 2) return; // 基底が1個以下なら何もせずに終了する

    Matrix B_star, U; // Gram_Schmidt の受け皿
    Gram_Schmidt(B, B_star, U);

    Vector B_norm = B_star.rowwise().squaredNorm(); // // ステップ 2: B_i = ||b*_i||^2 の計算
    int k_idx = 1; // 0-based index に注意し置き換えた
    while(k_idx < n){ 
        
        // ステップ 6: Size-reduce(k, j) (Size-reduce_partial では B と U のみが更新)
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);
        }
        
        // GSOベクトルの再計算と B_norm の更新
        // b*_k の再計算: b*_k = b_k - sum_{j=1}^{k-1} mu_{k,j} * b*_j
        Vector b_k_star = B.row(k_idx); 
        for (int j = 0; j < k_idx; j++) {
            // U(k_idx, j) * B_star.row(j) の減算
            Scalar mu_ld = U(k_idx, j);
            b_k_star -= mu_ld * B_star.row(j); 
        }
        B_star.row(k_idx) = b_k_star;
        B_norm(k_idx) = b_k_star.squaredNorm();

        // Lovász 条件チェック: B_k >= (delta - mu_{k,k-1}^2) * B_{k-1}
        if(B_norm(k_idx) >= ((delta - U(k_idx, k_idx - 1) * U(k_idx, k_idx - 1)) * B_norm(k_idx - 1))){
            // ステップ 10: Lovász条件を満たす場合
            k_idx++;
        }else{
            B.row(k_idx - 1).swap(B.row(k_idx)); // 基底ベクトルの交換: b_{k-1} <-> b_k
            //// ステップ 12: GSO情報の更新 (k は 1-based)
            GSOUpdate_LLL_partial(U, B_norm, k_idx + 1);
            k_idx = std::max(k_idx - 1, 1);
        }
    }
}