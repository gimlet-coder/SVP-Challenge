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

#include "Size-reduce_partial.hpp"


void Size_reduce_partial(Matrix &B, Matrix &U, const int i,const int j){
    // red_B と red_U を出力したい → B と U を更新して入出力両方に用いる

    //i, j が 1 <= j < i <= n を満たしているか確認
    if (i <= j || i < 0 || j < 0 || i >= B.rows() || j >= B.rows()) {
        throw std::out_of_range("Size_reduce: 不正なインデックス i または j です.");
    }
    // U(j,l) にアクセスするため U.rows() もチェックする
    if (i >= U.rows() || j >= U.cols() || j >= U.rows()) { 
         throw std::out_of_range("Size_reduce: U行列に対するインデックスが範囲外です。");
    }

    // --- ここまで引数が条件を満たしているかの確認 ---

    if(std::fabs(U(i, j)) > 0.5){
        double mu_ij = U(i, j);
        double q_double = std::round(mu_ij); // GSO係数の更新をするために double で q を計算する
        long long q_int = static_cast<long long>(q_double); // 基底ベクトルの更新用に long long (整数型) で q を計算
        B.row(i) -= q_int*B.row(j);

        for (int l = 0; l < j; l++){
            U(i, l) -= q_double*U(j, l);
        }
        U(i, j) -= q_double;
    }
}

#if 0

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    
    Matrix B(3, 3);
    B << 5, -3, -7,
         2, -7, -7,
         3, -10, 0;

    Matrix B_star;
    Matrix U;

    // それぞれ, Gram_Schmidt と Size-reduce の受け皿を用意

    try{ // エラーが発生しうる箇所を try で計算する
        Gram_Schmidt(B, B_star, U); // GSOを計算

        std::cout << "元の基底行列 B:\n" << B << std::endl << std::endl;
        std::cout << "GSO係数行列 U:\n" << U << std::endl << std::endl;

        int i_idx = 1, j_idx = 0;
        std::cout << "--- Size_reduce(i=" << i_idx << ", j=" << j_idx << ") を呼び出し ---" << std::endl;
        Size_reduce_partial(B, U, i_idx, j_idx);
        std::cout << "簡約後の基底行列 red_B:\n" << B << std::endl << std::endl;
        std::cout << "簡約後の係数行列 red_U:\n" << U << std::endl << std::endl;

    }catch (const std::out_of_range& oor) { // 範囲外エラー (out_of_range) が発生した場合
        std::cerr << "エラーが発生しました (範囲外アクセス): " << oor.what() << std::endl;
        return 1; // エラー終了
    } catch (const std::exception& e) { // その他の標準的なエラーが発生した場合
        std::cerr << "予期せぬエラーが発生しました: " << e.what() << std::endl;
        return 1; // エラー終了
    }
    /* memo:
        std::cerr エラーメッセージを表示する専門　(std::cout は通常のプログラムの出力)
        orr.what(), e.what() それぞれが表示するエラーメッセージの内容
    
    */

    return 0;
}

#endif

