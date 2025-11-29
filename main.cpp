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
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている



#include "main.hpp"


#if 0
// DeepLLLの動作確認用
int main() {
    SetConsoleOutputCP(65001); // このコードの出力文字コードを UTF-8 に強制
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    std::random_device rd; // 乱数のシード値を取得
    std::mt19937 gen(rd()); //MT法に基づいた乱数エンジン

    int min_val = 1'000'000, max_val = 10'000'000; // 乱数範囲

     int n = 4; // 格子次元
    
    std::uniform_int_distribution<int> dist(min_val, max_val); // 上記の乱数範囲に応じたランダムな整数を生成する関数
   
    /* Eigen::MatrixXi B_int(n, n);
    B_int.setZero();
    for (int i = 0; i < n; i++){
        B_int(i, i) = 1; // 対角成分を 1 に設定
        B_int(i, 0) = dist(gen); // 1列目に大きな乱数を設定
    }
   */
    Eigen::MatrixXi B_int(n, n);

    B_int << 388, 417, 417, -86,
            -672, -73, -121, 944,
            -689, 379, 724, 653,
            -179, 96, -24, 978;
    
   

    std::cout << "入力した基底行列 B:\n" << B_int << std::endl << std::endl;

    Scalar Norm_input = B_int.row(0).cast<Scalar>().squaredNorm();

    double delta = 0.75;



    Matrix B_long_double(n, n); 
    B_long_double = B_int.cast<typename Matrix::Scalar>(); // Matrix::Scalarはboost::multiprecision型

    try{
        LLL(B_long_double, delta);

        auto B_tmp_ld = B_long_double.cast<long double>();
        B_int = B_tmp_ld.array().round().cast<int>();
        
        std::cout << "簡約パラメータ delta = " << delta << " によってLLL基底簡約した B:\n" 
                    << std::fixed << std::setprecision(0) // 整数表示 (浮動小数点数演算を使用しているため)
                     << B_long_double << std::endl << std::endl;
    }
    catch(const std::out_of_range &orr){
        std::cerr << orr.what() << std::endl;
        return 1; // エラー終了
    }

    std::cout << "--------------------------------" << std::endl;
    std::cout << "簡約前後の最初のベクトルの長さ比較 (2乗ノルム):" << Norm_input <<std::endl;


// 出力(B_long_double)の最初の行の長さ
    std::cout << "LLL後: " << B_long_double.row(0).squaredNorm() << std::endl;
    return 0;
}

#endif

#if 1
// MLLLの動作確認用
int main() {
    SetConsoleOutputCP(65001); // このコードの出力文字コードを UTF-8 に強制
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    std::random_device rd; // 乱数のシード値を取得
    std::mt19937 gen(rd()); //MT法に基づいた乱数エンジン

    int min_val = 1'000'000, max_val = 10'000'000; // 乱数範囲

    // int n = 40; // 格子次元
    
    std::uniform_int_distribution<int> dist(min_val, max_val); // 上記の乱数範囲に応じたランダムな整数を生成する関数
    
    Eigen::MatrixXi B_int(5, 4);

    // 例5.2.1
    B_int << 388, 417, 417, -86,
            -672, -73, -121, 944,
            -689, 379, 724, 653,
            -179, 96, -24, 978,
            508, -705, 173, -343;
    
    std::cout << "入力した基底行列 B:\n" << B_int << std::endl << std::endl;

    Scalar Norm_input = B_int.row(0).cast<Scalar>().squaredNorm();

    double delta = 0.75;



    Matrix B_long_double(5, 4); 
    B_long_double = B_int.cast<typename Matrix::Scalar>(); // Matrix::Scalarはboost::multiprecision型

    try{
        MLLL(B_long_double, delta);

        auto B_tmp_ld = B_long_double.cast<long double>();
        B_int = B_tmp_ld.array().round().cast<int>();
        
        std::cout << "簡約パラメータ delta = " << delta << " によってMLLL基底簡約した B:\n" 
                    << std::fixed << std::setprecision(0) // 整数表示 (浮動小数点数演算を使用しているため)
                     << B_long_double << std::endl << std::endl;
    }
    catch(const std::out_of_range &orr){
        std::cerr << orr.what() << std::endl;
        return 1; // エラー終了
    }

    std::cout << "--------------------------------" << std::endl;
    std::cout << "簡約前後の最初のベクトルの長さ比較 (2乗ノルム):" << Norm_input <<std::endl;


// 出力(B_long_double)の最初の行の長さ
    std::cout << "MLLL後: " << B_long_double.row(0).squaredNorm() << std::endl;
    return 0;
}

#endif

