#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> 
#include <stdexcept> 

#include "lattice_types.hpp"

#include "DeepLLL.hpp"

// Size-reduce_partial 関数 (外部定義)
void Size_reduce_partial(Matrix &B, Matrix &U, const int i_in, const int j_in);
// Gram_Schmidt 関数 (外部定義)
void Gram_Schmidt(const Matrix &B, Matrix &B_star, Matrix &U);
// DeepLLL基底簡約で隣接するベクトルの交換後のGSO情報を更新する関数 (外部定義)
void GSOUpdate_DeepLLL_partial(Matrix &U, Vector &B_norm, const int i_in, const int k_in);



// アルゴリズム 7: DeepLLL 基底簡約アルゴリズム
void DeepLLL(Matrix &B, const double delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size_reduce: 不正なパラメータ δ です");
    }

    int n = B.rows();
    if(n < 2) return; 

    Matrix B_star, U; 
    Gram_Schmidt(B, B_star, U);

    Vector B_norm = B_star.rowwise().squaredNorm(); // B_i = ||b*_i||^2

    int k_idx = 1;
    long double delta_ld = static_cast<long double>(delta);

    // ステップ 4: while k <= n
    while(k_idx < n){
        
        // Step 6: Size-reduce(k, j)
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);
        }
        
        // Size-reduce後の B*_k と B_k の再計算 (B_norm(k_idx) の更新)
        Vector b_k_star = B.row(k_idx); 
        for (int j = 0; j < k_idx; j++) {
            long double mu_ld = U(k_idx, j).convert_to<long double>();
            b_k_star -= static_cast<Scalar>(mu_ld) * B_star.row(j);
        }
        B_star.row(k_idx) = b_k_star;
        B_norm(k_idx) = b_k_star.squaredNorm();
        
        Scalar C_scalar = B_norm(k_idx);
        int i_idx = 0;
        bool FLAG_DEEP_INSERTION = false;

        while (i_idx < k_idx){
            
            // C = ||pi_i(b_k)||^2 の効率的な計算 (Step 11)
            if(i_idx > 0) {
                 C_scalar -= U(k_idx, i_idx - 1) * U(k_idx, i_idx - 1) * B_norm(i_idx - 1);
            }
            
            // Deep Insertion 条件チェック: C < delta * B_i (Step 12)
            if(C_scalar.convert_to<long double>() >= delta_ld * B_norm(i_idx).convert_to<long double>()){ 
                i_idx++; 
            }
            else {
                // --- Deep Insertion 実行: B の更新 (Step 15-18) ---
                Vector temp_B = B.row(k_idx); 
                for (int j = k_idx; j >= i_idx + 1; j--){
                    B.row(j) = B.row(j - 1); // B の行をシフト
                }
                B.row(i_idx) = temp_B; // B の i_idx に旧 b_k を挿入

                // ★ GSO Update に備えた B_star の行交換とシフト (重要)
                // U と B_norm の更新には B_star は不要だが、次のループのために整合性を取る
                B_star.row(k_idx).swap(B_star.row(i_idx)); 
                for (int j = k_idx; j >= i_idx + 2; j--){
                    B_star.row(j).swap(B_star.row(j - 1)); // B_star の行をシフト
                }
                
                // ★ Step 19: GSO情報の更新 (U と B_norm のみ)
                GSOUpdate_DeepLLL_partial(U, B_norm, i_idx + 1, k_idx + 1); 
                
                // ★ B_star の再構築 (部分的な Gram-Schmidt)
                // 更新された U を使って、影響を受けた i_idx 行目以降の B_star の中身を再計算
                for (int i_update = i_idx; i_update < n; i_update++) {
                    B_star.row(i_update) = B.row(i_update); // b_i で初期化
                    for (int j_update = 0; j_update < i_update; j_update++) {
                        // 更新された U(i, j) と既存の B_star.row(j) を使用
                        long double mu_ld = U(i_update, j_update).convert_to<long double>();
                        B_star.row(i_update) -= static_cast<Scalar>(mu_ld) * B_star.row(j_update);
                    }
                }
                
                // k_idx の調整
                k_idx = std::max(i_idx + 1, 2) - 1; 
                FLAG_DEEP_INSERTION = true; 
                break;
            }
        }
        
        if(!FLAG_DEEP_INSERTION) k_idx++;
    }
}
