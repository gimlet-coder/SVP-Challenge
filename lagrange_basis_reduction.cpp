#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::round を使うために必要

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている

// アルゴリズム 2: Lagrange 基底簡約アルゴリズム
Matrix Lag_basis_red(const Matrix& L){

    Matrix Lag_L = L; 
    
    // --- ステップ 1-3: ||b_1|| <= ||b_2|| を満たすように交換 ---
    if(Lag_L.rows() < 2) return Lag_L; // 2次元未満ならそのまま返す
    
    // ||b_1|| > ||b_2|| の場合、交換
    if(Lag_L.row(0).squaredNorm() > Lag_L.row(1).squaredNorm()){ 
        Lag_L.row(0).swap(Lag_L.row(1)); // b_1 と b_2 を交換
    }
    
    // --- ステップ 4-7: do-while ループ ---
    do{
        // ステップ 5: q = -[<b_1, b_2> / ||b_1||^2] を計算
        Scalar frac = -Lag_L.row(0).dot(Lag_L.row(1)) / Lag_L.row(0).squaredNorm();
        int q = static_cast<int>(std::round(frac.template convert_to<long double>()));
        
        // v <- b_2 + q * b_1 を計算
        Vector v = Lag_L.row(1) + (static_cast<typename Matrix::Scalar>(q)) * Lag_L.row(0);
        
        // ステップ 6: 基底の取り直し (b_2 <- b_1, b_1 <- v)
        // Lag_L.row(1) には元の b_1、Lag_L.row(0) には新しい v が入る
        Lag_L.row(1).swap(Lag_L.row(0)); // b_1 と b_2 の位置を入れ替える
        Lag_L.row(0) = v; // 新しい b_1 (より短い v) を設定

    }while(Lag_L.row(0).squaredNorm() < Lag_L.row(1).squaredNorm()); // ||b_1|| < ||b_2|| の間繰り返す
    
    // --- ステップ 8: 最終的なソート: ||b_1|| <= ||b_2|| を満たすように交換 ---
    // do-while 終了時は ||b_1|| > ||b_2|| となっているので、再度交換して ||b_1|| <= ||b_2|| にする
    Lag_L.row(0).swap(Lag_L.row(1));
    return Lag_L;
}
