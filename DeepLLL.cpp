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
void DeepLLL(Matrix &B, const Scalar delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size_reduce: 不正なパラメータ δ です");
    }

    int n = B.rows();
    if(n < 2) return; 

    Matrix B_star, U; 
    Gram_Schmidt(B, B_star, U);

    Vector B_norm = B_star.rowwise().squaredNorm(); // B_i = ||b*_i||^2 step.2

    int k_idx = 1;

    // step.4: while k <= n
    while(k_idx < n){
        
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);// Step 6: Size-reduce(k, j)
        }
        
        // Size-reduce後の B*_k と B_k の再計算 (B_norm(k_idx) の更新)
        Vector b_k_star = B.row(k_idx);
        for (int j = 0; j < k_idx; j++) {
            b_k_star -= U(k_idx, j) * B_star.row(j);
        }
        B_star.row(k_idx) = b_k_star;
        B_norm(k_idx) = b_k_star.squaredNorm();
        
        Scalar C_scalar = B.row(k_idx).squaredNorm();
        int i_idx = 0;

        while(i_idx < k_idx){
            if(C_scalar >= delta * B_norm(i_idx)){
                C_scalar -= U(k_idx, i_idx) * U(k_idx, i_idx) * B_norm(i_idx);
                i_idx++;
            }else{
                Vector temp = B.row(k_idx);
                for (int j = k_idx; j > i_idx; j--){
                    B.row(j) = B.row(j - 1);
                }
                B.row(i_idx) = temp;
                GSOUpdate_DeepLLL_partial(U, B_norm, i_idx + 1, k_idx + 1);
                Gram_Schmidt(B, B_star, U);
                B_norm = B_star.rowwise().squaredNorm();
                k_idx = std::max(i_idx, 1) - 1;
            }
        }
        k_idx++;
    }
}