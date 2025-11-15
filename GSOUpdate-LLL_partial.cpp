#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar (=cpp_int) 型になっている



void GSOUpdate_LLL_partial(Matrix &U,Vector &B_norm, const int k){
    // 2 <= k <= n ( = U.rows()) を満たしているか確認
    if(k < 2 || U.rows() < k){
        throw std::out_of_range("Size_reduce: 不正なインデックス k です");
    }
    // --- 引数チェック終了 ---

    // 0-based に注意
    int k_idx = k - 1; // アルゴリズム6の表記と比較しやすいよう 1-based にした
    const int km1_idx = k_idx - 1; // k-1 の 0-based index

    // --- ステップ 1: nu, delta の計算 ---
    Scalar nu = U(k_idx, k_idx - 1); 
    Scalar delta = B_norm(k_idx) + nu * nu * B_norm(k_idx - 1);
    
    // --- ステップ 2: GSOノルム B_k, B_{k-1} および $\mu_{k, k-1}$ の更新
    long double nu_ld = nu.convert_to<long double>();
    long double B_norm_km1_ld = B_norm(km1_idx).convert_to<long double>();
    long double B_norm_k_ld = B_norm(k_idx).convert_to<long double>();
    long double delta_ld = delta.convert_to<long double>();

    // mu_{k, k-1} <- nu B_{k-1} / delta$
    U(k_idx, km1_idx) = static_cast<Scalar>(nu_ld * B_norm_km1_ld / delta_ld);

    // B_k <- B_k B_{k-1} / delta
    B_norm(k_idx) = static_cast<Scalar>(B_norm_k_ld * B_norm_km1_ld / delta_ld);

    // B_{k-1} <- \delta
    B_norm(km1_idx) = delta;

    // --- ステップ 3-5: mu_{k-1, j} と mu_{k, j} の交換 (j < k-1) (補題 2.3.7(2)) ---
    for (int j = 0; j < k - 2; j++){
        U(km1_idx, j).swap(U(k_idx, j)); // mu_{k-1, j} <-> mu_{k, j}
    }

    // --- ステップ 6-9: mu_{i, k-1} と mu_{i, k} の更新 (i > k) (補題 2.3.7(3)) ---
    for (int i = k; i < U.rows(); i++){

        Scalar temp = U(i, k_idx);
        long double U_ik_ld = U(i, k_idx).convert_to<long double>();
        long double U_ikm1_ld = U(i, km1_idx).convert_to<long double>();
        long double nu_Ukm1k_ld = U(k_idx, km1_idx).convert_to<long double>();

        long double result1_ld = U_ikm1_ld - nu_ld * U_ik_ld;
        U(i, k_idx) = static_cast<Scalar>(result1_ld);

        // ここでの U(i, k) は上記で更新された U(i, k_idx) の値を使用する
        long double U_ik_new_ld = U(i, k_idx).convert_to<long double>();
        long double result2_ld = temp.convert_to<long double>() + nu_Ukm1k_ld * U_ik_new_ld;
        U(i, km1_idx) = static_cast<Scalar>(result2_ld);
    }
}
