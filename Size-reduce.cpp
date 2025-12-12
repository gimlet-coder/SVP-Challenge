#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている

#include "Size-reduce.hpp"

// アルゴリズム 4: Size-reduce (完全なサイズ簡約)
void Size_reduce(IntMatrix &B, RealMatrix &U){
    const int n = B.rows(); // 行列 B の行数を取得して, ベクトルの本数 n として扱う
    if (n < 2){ // 基底ベクトルが1本以下の場合
        if(n > 0){
            RealMatrix B_star_temp;
            Gram_Schmidt(B, B_star_temp, U); // GSOを計算して U を返す
        }else{
            U.resize(0, 0); // 空行列の場合は U も空にする
        }
        return; // サイズ簡約の必要がないので終了
    }

    RealMatrix B_star; // GSOベクトル及びGSO係数の受け皿
    RealMatrix loop_U; // ループ内で用いる一時的なGSO係数行列
    Gram_Schmidt(B, B_star, loop_U);
    
    // ステップ 2-6: すべての 2 <= i <= n, 1 <= j < i に対して Size-reduce_partial を適用
    for (int i = 1; i < n; i++){
        for (int j = i - 1; j >= 0; j--){
            Size_reduce_partial(B, loop_U, i, j);
        }
    }
    // Size_reduce_partial では GSOベクトル B_star は変化しないが、最終的な U の値は更新された B に対して再計算する必要がある
    // (B_star が必要ないなら、Gram_Schmidt を呼び出すことで U を更新するだけでOK)
    Gram_Schmidt(B, B_star, U);
}