#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

// Eigenの型を使いやすく名前を変更した
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

// Eigen::MatrixXd は内部的には double 型

/*memo:
B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す
hoge.dot(huga) は hoge と huga の内積を表す
hoge.squaredNorm() はベクトルの長さの2乗を表す
*/

#include "LLL.hpp"

void LLL(Matrix &B, const double delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size-reduce 不正なインデックス delta です (1/4 < delta < 1)");
    }
    // --- 引数チェック終了 ---

    int n = B.rows();
    if(n < 2) return; // 基底が1個以下なら何もせずに終了する



    Matrix B_star, U; // Gram_Schmidt の受け皿
    Gram_Schmidt(B, B_star, U);

    Vector B_norm = B_star.rowwise().squaredNorm(); // hoge.rowwise() は hoge の行ごとの処理を一度に行う
    int k_idx = 1; // 0-based index に注意し置き換えた
    while(k_idx < n){
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);
        }
        if(B_norm(k_idx) >= ((delta - U(k_idx, k_idx - 1) * U(k_idx, k_idx - 1)) * B_norm(k_idx - 1))){
            k_idx++; //Lovasz条件を満たす場合
        }else{
            B.row(k_idx - 1).swap(B.row(k_idx)); // 基底ベクトルの交換 hoge(a).swap(hoge(b)) で hoge の a と b を交換する
            GSOUpdate_LLL_partial(U, B_norm, k_idx + 1);
            k_idx = std::max(k_idx - 1, 1);
        }
    }
}