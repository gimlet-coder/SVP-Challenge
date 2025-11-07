#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

// Eigenの型を使いやすく名前を変更した
using Vector = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;

// Eigen::MatrixXd は内部的には double 型

/*memo:
B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す
hoge.dot(huga) は hoge と huga の内積を表す
hoge.squaredNorm() はベクトルの長さの2乗を表す
*/

#include "DeepLLL.hpp"

void DeepLLL(Matrix &B, const double delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size_reduce: 不正なパラメータ δ です");
    }
    /* --- 引数チェック完了 --- */

    int n = B.rows();
    if(n < 2) return; // 基底が1個以下なら何もせずに終了する

    Matrix B_star, U; // Gram_Schmidt の受け皿
    Gram_Schmidt(B, B_star, U);

    Vector B_norm = B_star.rowwise().squaredNorm(); // hoge.rowwise() は hoge の行ごとの処理を一度に行う

    int k_idx = 1;

    

    while(k_idx < n){
        bool FLAG = false;
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);
        }
        double C = B_norm(k_idx);
        C = B.row(k_idx).squaredNorm();
        int i_idx = 0;
        while (i_idx < k_idx){
            C -= U(k_idx, i_idx) * U(k_idx, i_idx) * B_norm(i_idx);
            if(C >= delta * B_norm(i_idx)){
                i_idx++;
            }else{
                Vector temp = B.row(k_idx);
                for (int j = k_idx; j >= i_idx + 1; j--){
                    B.row(j) = B.row(j - 1);                    
                }
                B.row(i_idx) = temp;
                GSOUpdate_DeepLLL_partial(U, B_norm, i_idx + 1, k_idx + 1); // この関数の i, k に当たるものは 1-index であることに注意する
                k_idx = std::max(i_idx, 1);
                FLAG = true;
                break;
            }
        }
        if(!FLAG){
            k_idx++;
        }
    }

}