#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include <random> // 10^7 以上の疑似乱数を生成するために必要

#include <windows.h> //　文字化け対策
#include <iomanip> // 表示桁数を確保するのに使う std::fixed << std::setprecision(0)

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は long double 型になっている



#include "main.hpp"

// DeepLLL 関数 (外部定義)
void DeepLLL(Matrix &B, const double delta);

int main() {
    SetConsoleOutputCP(65001); // このコードの出力文字コードを UTF-8 に強制
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    std::random_device rd; // 乱数のシード値を取得
    std::mt19937 gen(rd()); //MT法に基づいた乱数エンジン

    int min_val = 10'000, max_val = 100'000; // 乱数範囲

    int n = 20; // 格子次元
    
    std::uniform_int_distribution<int> dist(min_val, max_val); // 上記の乱数範囲に応じたランダムな整数を生成する関数
    
    Eigen::MatrixXi B_int(n, n);

    B_int.setZero();
    for (int i = 0; i < n; i++){
        B_int(i, i) = 1; // 対角成分を 1 に設定
        B_int(i, 0) = dist(gen); // 1列目に大きな乱数を設定
    }
    
    std::cout << "入力した基底行列 B:\n" << B_int << std::endl << std::endl;
    double delta = 0.99;



    Matrix B_long_double(n, n); 
    B_long_double = B_int.cast<typename Matrix::Scalar>(); // Matrix::Scalarはboost::multiprecision型

    try{
        DeepLLL(B_long_double, delta);

        Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> B_tmp_ld = B_long_double.cast<long double>();
        B_int = B_long_double.cast<long double>().array().round().cast<int>();
        
        std::cout << "簡約パラメータ delta = " << delta << " によってDeepLLL基底簡約した B:\n" 
                    << std::fixed << std::setprecision(0) // 整数表示 (浮動小数点数演算を使用しているため)
                     << B_long_double << std::endl << std::endl;
    }
    catch(const std::out_of_range &orr){
        std::cerr << orr.what() << std::endl;
        return 1; // エラー終了
    }
    return 0;
}