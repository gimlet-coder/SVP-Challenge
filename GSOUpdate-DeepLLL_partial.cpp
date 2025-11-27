#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> 
#include <stdexcept> 

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている
void GSOUpdate_DeepLLL_partial(Matrix &U, Vector &B_norm, const int i_in, const int k_in){
    // i_in, k_in は 1-based index (教科書の表記に合わせる)
    int n = U.rows();
    if(i_in < 1 || k_in <= i_in || n < k_in ){
        throw std::out_of_range("GSOUpdate_DeepLLL_partial: 不正なインデックス i または k です");
    }
    const int i = i_in - 1, k = k_in - 1; // 0-based index に変換
    
    // P, D ベクトルの計算に必要な初期化
    Vector P(n), D(n); 

    // --- ステップ 1-5: P_j と D_j (||pi_j(b_k)||^2) の計算 (補題 2.4.2) ---
    D(k) = B_norm(k); // D_k = ||pi_k(b_k)||^2 = ||b_k^*||^2 = B_k
    P(k) = B_norm(k);
    // j = k - 1 から i まで、D_j = D_{j+1} + mu_{k,j}^2 * B_j を計算
    for (int j = k - 1; j >= i; j--){
        P(j) = U(k, j) * B_norm(j); // P_j = mu_k,j * B_j (アルゴリズム8のP_jとは定義が異なるが、ここでは mu_{k,j} * B_j を表す)
        D(j) = D(j + 1) + U(k, j) * P(j); // D_j = ||pi_j(b_k)||^2
    }
    
    Vector S(n); // S_l = sum_{s=j}^{k} nu_s * mu_{l,s} * B_s (補題 2.4.3(1)の計算用)
    S.setZero(); 

    // --- ステップ 13-14: U(l, j) の更新 (j > i) (補題 2.4.3(1) の後半) ---
    for (int j = k; j >= i + 1; j--){
        // T = mu_{k,j-1} / D_j。 Deep Insertion 前は mu_{k,j-1} は存在しないので 0
        Scalar T;
        if(j - 1 >= 0) {
            T = U(k, j - 1) / D(j);
        }
        
        
        // l > k の更新 (l = n - 1 から k)
        for (int l = n - 1; l >= k + 1; l--){
            
            // U(l, j) の更新: nu_{l,j} <- mu_{l,j-1} - T * S_l
            S(l) += U(l, j) * P(j);
            U(l, j) = U(l, j - 1) - T * S(l);
        }
        
        // k >= l > j の更新 (l = k から j + 1)
        for (int l = k; l >= j + 1; l--){

            Scalar mu_val;

            if(l - 1 == j){
                mu_val = 1;
            }else{
                mu_val = U(l - 1, j);
            }

            S(l) += mu_val * P(j);
            U(l, j) = U(l - 1, j - 1) - T * S(l);
        }
    }
    
    // --- ステップ 16-23: mu_l,i の更新 (補題 2.4.3(2)) ---
    
    Scalar T = 1.0 / D(i);

    // l > k の更新 (l = n - 1 から k + 1)
    for (int l = n - 1; l >= k + 1; l--){
        U(l, i) = T * (S(l) + U(l, i) * P(i));
    }
    
    // k >= l > i の更新 (l = k - 1 から i + 2)
    for (int l = k; l >= i + 2; l--){
        U(l, i) = T * (S(l) + U(l - 1, i) * P(i));
    }

    U(i + 1, i) = T * P(i);

    // --- ステップ 24-30: mu_l,j のシフト (j < i) (補題 2.4.3(3), (4)) ---
    for (int j = 0; j < i; j++){
        Scalar temp = U(k, j);
        for (int l = k; l >= i + 1; l--){
            U(l, j) = U(l - 1, j);
        }
        U(i, j) = temp;
    }
    
    // --- ステップ 31-34: B_norm の更新 (補題 2.4.2) ---
    for (int j = k; j >= i + 1; j--){
        // B_j <- D_j * B_{j-1} / D_{j-1}
        B_norm(j) = D(j) * B_norm(j - 1) / D(j - 1);
    }
    B_norm(i) = D(i); // B_i <- D_i = ||pi_i(b_k)||^2
}