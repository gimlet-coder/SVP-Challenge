#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::round を使うために必要

#include "lattice_types.hpp"
// Eigen::MatrixXd は内部的には double 型

// 2次元格子 L の基底{b_1, b_2} を Lagrange 基底簡約して返す

Matrix Lag_basis_red(const Matrix& L){

    Matrix Lag_L = L; //Lag_L を L で初期化する 行数がそれぞれ一致しているためシンプルにこれでOK

    /*memo:
        B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
        B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す

        hoge.dot(huga) は hoge と huga の内積を表す
        hoge.squaredNorm() はベクトルの長さの2乗を表す
    */

    if(Lag_L.row(0).squaredNorm() > Lag_L.row(1).squaredNorm()){
        Vector v = Lag_L.row(0);
        Lag_L.row(0) = Lag_L.row(1); 
        Lag_L.row(1) = v;
    }
    do{
        double fraq = -Lag_L.row(0).dot(Lag_L.row(1)) / Lag_L.row(0).squaredNorm();
        int q =  std::round(fraq);
        Vector v = Lag_L.row(1) + q * Lag_L.row(0);
        Lag_L.row(1) = Lag_L.row(0); Lag_L.row(0) = v;
    }while(Lag_L.row(0).squaredNorm() < Lag_L.row(1).squaredNorm());
    Vector temp_v = Lag_L.row(0); Lag_L.row(0) = Lag_L.row(1); Lag_L.row(1) = temp_v;
    return Lag_L;
}

#if 0
int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    
    Matrix L(2, 4);
    L << 230, -651, 609, -366,
         301, -852, 797, -479;

    Matrix Lag_L = Lag_basis_red(L);

    std::cout << "imput Lattice L \n" << L << std::endl;
    std::cout << "Lagrangian reduced basis for lattice L \n" << Lag_L << std::endl;
    

    return 0;
}

#endif