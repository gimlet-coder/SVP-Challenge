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


// アルゴリズム 5: LLL 基底簡約アルゴリズム
void LLL(IntMatrix &B, const Real delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size-reduce 不正なインデックス delta です (1/4 < delta < 1)");
    }
    // --- 引数チェック終了 ---

    int n = B.rows();
    if(n < 2) return; // 基底が1個以下なら何もせずに終了する

    RealMatrix B_star, U; // Gram_Schmidt の受け皿
    Gram_Schmidt(B, B_star, U);

    RealVector B_norm = B_star.rowwise().squaredNorm(); // // ステップ 2: B_i = ||b*_i||^2 の計算
    int k_idx = 1; // 0-based index に注意し置き換えた
    while(k_idx < n){ 
        
        // ステップ 6: Size-reduce(k, j) (Size-reduce_partial では B と U のみが更新)
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);
        }
        
        // Lovász 条件チェック: B_k >= (delta - mu_{k,k-1}^2) * B_{k-1}
        if(B_norm(k_idx) >= ((delta - U(k_idx, k_idx - 1) * U(k_idx, k_idx - 1)) * B_norm(k_idx - 1))){
            // ステップ 10: Lovász条件を満たす場合
            k_idx++;
        }else{
            B.row(k_idx - 1).swap(B.row(k_idx)); // 基底ベクトルの交換: b_{k-1} <-> b_k
            //// ステップ 12: GSO情報の更新 (k は 1-based)
            GSOUpdate_LLL_partial(U, B_norm, k_idx);
            k_idx = std::max(k_idx - 1, 1);
        }
    }
}