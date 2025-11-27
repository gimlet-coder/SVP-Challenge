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
    // 1 <= k < n ( = U.rows()) を満たしているか確認
    if(k < 1 || U.rows() <= k){
        throw std::out_of_range("Size_reduce: 不正なインデックス k です");
    }
    // --- 引数チェック終了 ---
    Scalar nu = U(k, k - 1);
    Scalar delta = B_norm(k) + nu * nu * B_norm(k - 1);
    U(k, k - 1) = nu * B_norm(k - 1) / delta;
    B_norm(k) *= B_norm(k - 1) / delta;
    B_norm(k - 1) = delta;

    for (int j = 0; j <= k - 2; j++){
        U(k - 1, j).swap(U(k, j));
    }
    for (int i = k + 1; i < U.rows(); i++){
        Scalar temp = U(i, k);
        U(i, k) = U(i, k - 1) - nu * temp;
        U(i, k - 1) = temp + U(k, k - 1) * U(i, k);
    }
    
}
