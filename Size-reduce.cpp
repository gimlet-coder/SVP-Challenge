#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include "lattice_types.hpp"

// Eigen::MatrixXd は内部的には double 型

/*memo:
B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す
hoge.dot(huga) は hoge と huga の内積を表す
hoge.squaredNorm() はベクトルの長さの2乗を表す
*/

#include "Size-reduce.hpp"

void Size_reduce(Matrix &B, Matrix &U){
    const int n = B.rows(); // 行列 B の行数を取得して, ベクトルの本数 n として扱う
    if (n < 2){ // 基底ベクトルが1本以下の場合
        if(n > 0){
            Matrix B_star_temp;
            Gram_Schmidt(B, B_star_temp, U); // GSOを計算して U を返す
        }else{
            U.resize(0, 0); // 空行列の場合は U も空にする
        }
        return; // サイズ簡約の必要がないので終了
    }

    Matrix B_star; // GSOベクトル及びGSO係数の受け皿
    Matrix loop_U; // ループ内で用いる一時的なGSO係数行列
    Gram_Schmidt(B, B_star, loop_U);
    
    for (int i = 1; i < n; i++){
        for (int j = i - 1; j >= 0; j--){
            Size_reduce_partial(B, loop_U, i, j);
        }
    }
    Gram_Schmidt(B, B_star, U);
}

#if 0
int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    Matrix B(3, 3);
    B << 5, -3, -7,
         2, -7, -7,
         3, -10, 0;

    std::cout << "元の基底行列 B:\n" << B << std::endl << std::endl;

    Matrix U_result;

    try{
        Size_reduce(B, U_result);

        std::cout << "--- サイズ簡約後の基底行列 B ---" << std::endl;
        std::cout << B << std::endl;

        std::cout << "\n--- 対応する GSO 係数行列 U ---" << std::endl;
        std::cout << U_result << std::endl;

    } catch (const std::out_of_range& oor) {
        std::cerr << "エラーが発生しました (範囲外アクセス): " << oor.what() << std::endl;
        return 1; // エラー終了
    } catch (const std::exception& e) {
        std::cerr << "予期せぬエラーが発生しました: " << e.what() << std::endl;
        return 1; // エラー終了
    }

    return 0;
}

#endif