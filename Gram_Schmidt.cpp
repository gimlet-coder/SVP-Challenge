#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている


void Gram_Schmidt(const IntMatrix& B, RealMatrix& B_star, RealMatrix& U){
    const int n = B.rows(); // 行列 B の行数を取得して, ベクトルの本数 n として扱う
    const int m = B.cols(); // 行列 B の列数を取得して, ベクトルの次元 m として扱う

    B_star.resize(n, m); U.resize(n, n); // B_star と U のサイズを指定
    U.setZero(); // U の要素を初期化



    for (int i = 0; i < n; i++){
            // B.row(i) は Integer なので、計算用に Real にキャストする
        RealVector b_i_real = B.row(i).cast<Real>();
        B_star.row(i) = b_i_real;
        
        for (int j = 0; j < i; j++){
            
            Real B_j_norm_sq = B_star.row(j).squaredNorm();
            if (abs(B_j_norm_sq) < MULTI_PRECISION_EPSILON){
                U(i, j) = 0; // 線形従属の場合
            } else {
                // --- ステップ 4: mu_i,j の計算 ---
                // mu_{i,j} = <b_i, b*_j> / ||b*_j||^2
                Real dot_product = b_i_real.dot(B_star.row(j));
                Real mu = dot_product / B_j_norm_sq;
                U(i, j) = mu;

                // ステップ 5: b*_i の更新
                B_star.row(i) -= mu * B_star.row(j);
            }
        }
        U(i, i) = 1; // GSO係数行列 U の対角成分は 1 となる
    }
}
