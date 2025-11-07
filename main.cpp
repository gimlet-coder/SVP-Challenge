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


#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は long double 型になっている


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

     Eigen::MatrixXi B_int(4, 4);

    B_int << 84,  3, 34, 17,
             20, 48, 66, 19,
             69, 14, 63, 78,
             28, 72, 36, 57;
    
    std::cout << "入力した基底行列 B:\n" << B_int << std::endl << std::endl;
    double delta = 0.75;

    Matrix B_long_double = B_int.cast<long double>(); 
    try{
        DeepLLL(B_long_double, delta);

        B_int = B_long_double.array().round().cast<int>(); // round で数値誤差を丸めて int にキャストして値を戻す

        std::cout << "簡約パラメータ delta = " << delta << " によってDeepLLL基底簡約した B:\n" << B_long_double << std::endl << std::endl;
    }
    catch(const std::out_of_range &orr){
        std::cerr << orr.what() << std::endl;
        return 1; // エラー終了
    }
    return 0;
}




#if 0
// LLL の動作確認
int main() {
    SetConsoleOutputCP(65001); // このコードの出力文字コードを UTF-8 に強制
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    std::random_device rd; // 乱数のシード値を取得
    std::mt19937 gen(rd()); //MT法に基づいた乱数エンジン

    int min_val = 10'000'000, max_val = 100'000'000; // 10^7 ~ 10^8 と指定

    std::uniform_int_distribution<int> dist(min_val, max_val); // 10^7~10^8 のランダムな整数を生成する関数
    
    Eigen::MatrixXi B_int(20, 20);

    B_int.setZero(); // B の初期化
    for (int i = 0; i < 20; i++){
        B_int(i, i) = 1;
        B_int(i, 0) = dist(gen);
    }
    
    std::cout << "入力した基底行列 B:\n" << B_int << std::endl << std::endl;
    double delta = 0.99;

    Eigen::MatrixXd B_double = B_int.cast<double>(); // LLL関数に入れるため double にキャストしたコピーを作成
    try{
        LLL(B_double, delta);
        
        B_int = B_double.array().round().cast<int>(); // round で数値誤差を丸めて int にキャストして値を戻す

        std::cout << "簡約パラメータ delta = " << delta << " によってLLL基底簡約した B:\n" << B_double << std::endl << std::endl;
    }
    catch(const std::out_of_range &orr){
        std::cerr << orr.what() << std::endl;
        return 1; // エラー終了
    }
    return 0;
}

#endif