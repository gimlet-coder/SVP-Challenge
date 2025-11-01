#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include <windows.h> //　文字化け対策

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


#include "main.hpp"

int main() {
    SetConsoleOutputCP(65001); // このコードの出力文字コードを UTF-8 に強制
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    
    Matrix B(3, 3);
    B << 9, 2, 7,
         8, 6, 1,
         3, 2, 6;
    
    std::cout << "入力した基底行列 B:\n" << B << std::endl << std::endl;
    double delta = 0.99;

    try{
        LLL(B, delta);
        std::cout << "簡約パラメータ delta = " << delta << " によってLLL基底簡約した B:\n" << B << std::endl << std::endl;
    }
    catch(const std::out_of_range &orr){
        std::cerr << orr.what() << std::endl;
        return 1; // エラー終了
    }

    return 0;
}