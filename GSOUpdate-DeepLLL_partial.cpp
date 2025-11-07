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

void GSOUpdate_DeepLLL_partial(Matrix &U, Vector &B_norm, const int i_in, const int k_in){
    int n = U.rows();
    /*--- i と k の引数チェック --- */
    if(i_in < 1 || k_in <= i_in || n < k_in ){
        throw std::out_of_range("Size_reduce: 不正なインデックス i または k です");
    }
    /* --- 引数チェック完了 --- */

    Vector P = B_norm, D = B_norm;
    
    const int i = i_in - 1, k = k_in - 1; // 0-index と 1-index の混同を防ぐためにずらした

    for (int j = k - 1; j >= i - 1; j--){
        P(j) = U(k, j) * B_norm(j);
        D(j) = D(j + 1) + U(k, j) * P(j);
    }
    Vector S(n);
    S.setZero(); // 初期化 アルゴリズムの6行目の処理

    for (int j = k - 1; j >= i; j--){
        double T = U(k, j - 1) / D(j);
        for (int l = n - 1; l >= k + 1; l--){
            S(l) += U(l, j) * P(j);
            U(l, j) = U(l, j - 1) - T * S(l);
        }
        for (int l = k - 1; l >= j + 1; l--){
            S(l) += U(l - 1, j) * P(j);
            U(l, j) = U(l - 1, j - 1) - T * S(l);
        }
    }
    double T = 1.0/D(i);
/* ここまででアルゴリズムの 16行目の処理 */
    for (int l = n - 1; l >= k + 1; l--){
        U(l, i) = T * (S(l) + U(l, i) * P(i - 1));
    }
    for (int l = k - 1; l >= i + 1; l--){
        U(l, i) = T * (S(l) + U(l - 1, i) * P(i - 1));
    }
    U(i, i - 1) = T * P(i - 1);
    for (int j = 0; j <= i - 2; j++){
        double temp = U(k, j);
        for (int l = k - 1; l >= i; l--){
            U(l + 1, j) = U(l, j);
        }
        U(i, j) = temp;
    }
    for (int j = k - 1; j >= i; j--){
        B_norm(j) = D(j) * B_norm(j - 1) / D(j - 1);
    }
    B_norm(i - 1) = D(i - 1);
}
