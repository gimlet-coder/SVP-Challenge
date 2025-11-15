#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include "lattice_types.hpp"

#include "Size-reduce_partial.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている

void Size_reduce_partial(Matrix &B, Matrix &U, const int i_in,const int j_in){
    // red_B と red_U を出力したい → B と U を更新して入出力両方に用いる

    //i, j が 1 <= j < i <= n を満たしているか確認
    if (i_in <= j_in || i_in < 0 || j_in < 0 || i_in >= B.rows() || j_in >= B.rows()) {
        throw std::out_of_range("Size_reduce: 不正なインデックス i または j です.");
    }
    // U(j,l) にアクセスするため U.rows() もチェックする
    if (i_in >= U.rows() || j_in >= U.cols() || j_in >= U.rows()) { 
         throw std::out_of_range("Size_reduce: U行列に対するインデックスが範囲外です。");
    }

    // --- ここまで引数が条件を満たしているかの確認 ---
    const int i = i_in, j = j_in; // 0-based index

    long double mu_ij_ld = U(i, j).convert_to<long double>();

    if(std::fabs(mu_ij_ld) > 0.5){ 
        
        // --- ステップ 2: q = [mu_{i,j}] を計算 ---
        long double q_ld = std::round(mu_ij_ld);
        // q は整数係数 (Scalar) で定義
        Scalar q = static_cast<Scalar>(q_ld); 
        
        // --- ステップ 2: 基底ベクトル更新: b_i <- b_i - q * b_j (Scalar (多倍長整数) 演算) ---
        B.row(i) -= q * B.row(j);
        
        // --- ステップ 4: GSO係数 mu_{i,l} の更新: mu_{i,l} <- mu_{i,l} - q * mu_{j,l} ---
        for (int l = 0; l < j; l++){
            // U の要素を long double で計算し、Scalar に戻す
            long double U_il_ld = U(i, l).convert_to<long double>();
            long double U_jl_ld = U(j, l).convert_to<long double>();
            
            long double result_ld = U_il_ld - q_ld * U_jl_ld;
            U(i, l) = static_cast<Scalar>(result_ld);
        }
        
        // --- U(i, j) の更新: mu_{i,j} <- mu_{i,j} - q ---
        // (この時点での U(i, j) の値は、更新された b_i に対応する新しい mu_{i,j} の値)
        long double U_ij_ld = U(i, j).convert_to<long double>();
        U_ij_ld -= q_ld; 
        U(i, j) = static_cast<Scalar>(U_ij_ld); 
    }
}