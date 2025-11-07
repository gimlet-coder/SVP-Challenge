#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> 
#include <stdexcept> 

// Eigenの型を long double に切り替える
using Vector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;

// DeepLLLのGSO更新 (アルゴリズム8のロジックに忠実な最終修正案)
void GSOUpdate_DeepLLL_partial(Matrix &U, Vector &B_norm, const int i_in, const int k_in){
    int n = U.rows();
    if(i_in < 1 || k_in <= i_in || n < k_in ){
        throw std::out_of_range("GSOUpdate_DeepLLL_partial: 不正なインデックス i または k です");
    }
    const int i = i_in - 1, k = k_in - 1; 
    
    // P, D ベクトルの計算に必要な初期化
    Vector P(n), D(n); 

    D(k) = B_norm(k);
    
    for (int j = k - 1; j >= i; j--){
        P(j) = U(k, j) * B_norm(j); // P_j = mu_k,j * B_j
        D(j) = D(j + 1) + U(k, j) * P(j); // D_j = ||pi_j(b_k)||^2
    }
    
    Vector S(n);
    S.setZero(); 

    // mu_l,j の更新 (j > i)
    for (int j = k; j >= i + 1; j--){
        // T = mu_k,j-1 / D_j。 Deep Insertion前は mu_k,j-1 は存在しないので 0
        long double T = 0.0L;
        if(j - 1 >= 0) {
            T = U(k, j - 1) / D(j);
        }
        
        // l > k の更新
        for (int l = n - 1; l >= k + 1; l--){
            S(l) += U(l, j) * P(j);
            U(l, j) = U(l, j - 1) - T * S(l);
        }
        
        // k >= l > j の更新
        for (int l = k; l >= j + 1; l--){
            S(l) += U(l - 1, j) * P(j);
            U(l, j) = U(l - 1, j - 1) - T * S(l);
        }
    }
    
    // mu_l,i の更新 (補題 2.4.3(2))
    long double T_div_D_i = 1.0 / D(i);

    // l > k の更新
    for (int l = n - 1; l >= k + 1; l--){
        U(l, i) = T_div_D_i * (S(l) + U(l, i) * P(i)); 
    }
    
    // k >= l > i の更新
    for (int l = k; l >= i + 1; l--){
        U(l, i) = T_div_D_i * (S(l) + U(l - 1, i) * P(i)); 
    }
    
    // mu_l,j のシフト (j < i)
    for (int j = 0; j <= i - 1; j++){
        double temp = U(k, j);
        for (int l = k; l >= i + 1; l--){
            U(l, j) = U(l - 1, j);
        }
        U(i, j) = temp;
    }
    
    // B_norm の更新 (補題 2.4.2)
    for (int j = k; j >= i + 1; j--){
        B_norm(j) = D(j) * B_norm(j - 1) / D(j - 1);
    }
    B_norm(i) = D(i);
}