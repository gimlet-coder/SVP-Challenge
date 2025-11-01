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

void GSOUpdate_LLL_partial(Matrix &U,Vector &B_norm, const int k){
    // 2 <= k <= n ( = U.rows()) を満たしているか確認
    if(k < 2 || U.rows() < k){
        throw std::out_of_range("Size_reduce: 不正なインデックス k です");
    }
    // --- 引数チェック終了 ---

    // 0-based に注意
    int k_idx = k - 1; // アルゴリズム6の表記と比較しやすいよう 1-based にした

    double nu = U(k_idx, k_idx - 1); 
    double delta = B_norm(k_idx) + nu * nu * B_norm(k_idx - 1);
    

    U(k_idx, k_idx - 1) = nu * B_norm(k_idx - 1)/delta;
    B_norm(k_idx) *=B_norm(k_idx - 1)/delta;
    B_norm(k_idx - 1) = delta;
    
    for (int j = 0; j < k - 2; j++){
        double temp = U(k_idx - 1, j);
        U(k_idx - 1, j) = U(k_idx, j);
        U(k_idx, j) = temp;
    }
    for (int i = k; i < U.rows(); i++){
        double temp = U(i, k_idx);
        U(i, k_idx) = U(i, k_idx - 1) - nu*temp;
        U(i, k_idx - 1)  = temp + U(k_idx, k_idx - 1) * U(i, k_idx);
    }
}
