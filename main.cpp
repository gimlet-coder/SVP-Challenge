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

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    
    return 0;
}