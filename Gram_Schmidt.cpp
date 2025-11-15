#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている


void Gram_Schmidt(const Matrix& B, Matrix& B_star, Matrix& U){
    const int n = B.rows(); // 行列 B の行数を取得して, ベクトルの本数 n として扱う
    const int m = B.cols(); // 行列 B の列数を取得して, ベクトルの次元 m として扱う

    B_star.resize(n, m); U.resize(n, n); // B_star と U のサイズを指定
    U.setZero(); // U の要素を初期化

    for (int i = 0; i < n; i++){
        B_star.row(i) = B.row(i); // b_i^* を b_i で初期化する
        
        for (int j = 0; j < i; j++){
            
            Scalar B_j_norm_scalar = B_star.row(j).squaredNorm();
            if (B_j_norm_scalar.compare(MULTI_PRECISION_EPSILON) < 0){
                U(i, j) = Scalar(0); // 線形従属の場合
            } else {
                // --- ステップ 4: mu_i,j の計算 ---
                // mu_{i,j} = <b_i, b*_j> / ||b*_j||^2
                Scalar dot_product = B.row(i).dot(B_star.row(j)); 
                Scalar mu = dot_product; 
                mu /= B_j_norm_scalar; 
                U(i, j) = mu;

                // ステップ 5: b*_i の更新
                B_star.row(i) -= mu * B_star.row(j);
            }
        }
        U(i, i) = Scalar(1); // GSO係数行列 U の対角成分は 1 となる
    }
}
