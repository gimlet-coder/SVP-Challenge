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

void Size_reduce_partial(Matrix &B, Matrix &U, const int i, const int j){
    // red_B と red_U を出力したい → B と U を更新して入出力両方に用いる

    //i, j が 1 <= j < i <= n を満たしているか確認
    if (i <= j || i < 0 || j < 0 || i >= B.rows() || j >= B.rows()) {
        throw std::out_of_range("Size_reduce: 不正なインデックス i または j です.");
    }
    // U(j,l) にアクセスするため U.rows() もチェックする
    if (i >= U.rows() || j >= U.cols() || j >= U.rows()) { 
         throw std::out_of_range("Size_reduce: U行列に対するインデックスが範囲外です。");
    }

    // --- ここまで引数が条件を満たしているかの確認 ---

    Scalar mu = U(i, j);
    if(abs(mu) >= Scalar(0.5)){
        Scalar q = round(mu);
        B.row(i) -= q * B.row(j);
        for (int l = 0; l <= j; l++){
            U(i, l) -= q * U(j, l); // GSO係数の更新
        }
    }
}
