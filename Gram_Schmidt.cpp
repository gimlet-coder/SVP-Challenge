#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は long double 型になっている

    /* const Matrix& B 入力用の行列だから中身も変えないしコピーも作らない */
    /* Matrix& B_star 出力用の行列 計算後のGSOベクトルを書き込むから const はつけない */
    /* Matrix& U 出力用の行列 GSOの係数を書き込むから const はつけない */

void Gram_Schmidt(const Matrix& B, Matrix& B_star, Matrix& U){
    const int n = B.rows(); // 行列 B の行数を取得して, ベクトルの本数 n として扱う
    const int m = B.cols(); // 行列 B の列数を取得して, ベクトルの次元 m として扱う

    B_star.resize(n, m); U.resize(n, n); // B_star と U のサイズを指定
    U.setZero(); // U の要素を初期化

    //ここまで準備 ここから下が本番
    /*memo:
        B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
        B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す

        hoge.dot(huga) は hoge と huga の内積を表す
        hoge.squaredNorm() はベクトルの長さの2乗を表す
    */

    for (int i = 0; i < n; i++){
        B_star.row(i) = B.row(i); //b_i^* を b_i で初期化する
        for (int j = 0; j < i; j++){
            long double mu = B_star.row(i).dot(B_star.row(j)) / B_star.row(j).squaredNorm(); // これが μ_(i, j)
            U(i, j) = mu; // 計算した μ_(i, j) を U に保存する
            B_star.row(i) -= mu * B_star.row(j);
        }
        U(i, i) = 1.0L; // GSO係数行列 U の対角成分は 1 となる
    }
}

#if 0
int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL); // 早く動かすおまじない

    Matrix B(3, 3);
    B << 5, -3, -7,
         2, -7, -7,
         3, -10, 0;

    Matrix B_star;
    Matrix U;

    Gram_Schmidt(B, B_star, U);

    std::cout << "input Matrix B:\n" << B << std::endl << std::endl;
    std::cout << "GSO Vectors B* (Orthogonal Basis):\n" << B_star << std::endl << std::endl;
    std::cout << "GSO Coefficients U (mu_ij):\n" << U << std::endl;

    return 0;
}
#endif