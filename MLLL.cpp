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

#include "MLLL.hpp"


void MLLL(Matrix &B, const Scalar delta){
    if(delta < 0.25 || 1 <= delta){
        throw std::out_of_range("Size-reduce 不正なインデックス delta です (1/4 < delta < 1)");
    }
    /* 引数チェック終了 */

    /* Step 1 */
    int h = B.rows(), z = h, g = 0;
    Matrix B_star, U; // GSOを求める際の受け皿
    Gram_Schmidt(B, B_star, U);

    Vector B_sq_norm(h); // 疑似コードでいうところの B にあたるもの
    for (int i = 0; i < h; i++){
        B_sq_norm(i) = B_star.row(i).squaredNorm();
    }

    while(g < z - 1){
        if(B_sq_norm(g) < 1e-20){ // ノルムの2乗が十分小さければ0ベクトルであると判定する(浮動小数点の誤差対策)
            if(g < z - 1){
                B.row(g).swap(B.row(z - 1));
            }
            z--;
            Gram_Schmidt(B, B_star, U);
            for (int i = 0; i < h; i++){
                B_sq_norm(i) = B_star.row(i).squaredNorm(); // GSO情報を更新するたびに B も更新することを忘れないように注意
            }
            continue;
        }
        
        if(g == 0){
            g = 1;
            continue;
        }
        int l = g, k = g;
        bool FLAG_startagain = false;
        while(k <= l && !FLAG_startagain){
            Size_reduce_partial(B, U, k, k - 1);
            Scalar nu = U(k, k - 1), B_proj = B_sq_norm(k) + nu*nu * B_sq_norm(k - 1);
            if(B_proj >= delta * B_sq_norm(k - 1)){
            /* step.15 */
                for(int j = k - 2; j >= 0; j--){
                    Size_reduce_partial(B, U, k, j);
                }
                k++;
            }else{
                /* step.17 */
                if(B_sq_norm(k) < 1e-20){ //上記同様浮動小数の誤差対策で判定基準を設ける
                    if(k < z - 1){
                        B.row(k).swap(B.row(k - 1));
                        Gram_Schmidt(B, B_star, U);
                        for(int i = 0; i < h; i++){
                            B_sq_norm(i) = B_star.row(i).squaredNorm();
                        }
                    }
                    z--;
                    g = k;
                    FLAG_startagain = true;
                }else{
                    B_star.row(k).swap(B_star.row(k - 1));
                    B.row(k).swap(B.row(k - 1));
                    for(int j = 1; j <= k - 2; j++){
                        std::swap(U(k, j - 1), (U(k - 1, j - 1)));
                    }
                    if(B_proj > 1e-20){ // 閾値でノンゼロチェックをする
                        /* step.24 */
                        if(B_sq_norm(k - 1) < 1e-20){
                            B_sq_norm(k - 1) = B_proj;
                            B_star.row(k - 1) = B_star.row(k);
                            B_star.row(k - 1) *= nu;
                            U(k, k - 1) = 1/nu;
                        }else{
                            /* step. 28 */
                            Scalar t = B_sq_norm(k - 1) / B_proj;
                            Scalar nu_old = nu;

                            U(k, k - 1) = nu_old * t;
                            Vector w = B_star.row(k - 1);

                            B_star.row(k - 1) = B_star.row(k) + nu_old * w.transpose();
                            B_sq_norm(k - 1) = B_proj;
                            if(k <= l){
                                B_star.row(k) = -nu_old * B_star.row(k) + (B_sq_norm(k) / B_proj) * w.transpose();
                                B_sq_norm(k) *= t;
                                B_sq_norm(k - 1) = B_proj;
                            }
                            for (int i = k + 1; i <= l; i++){
                                Scalar t_u = U(i, k);
                                U(i, k) = U(i, k - 1) - nu * t_u;
                                U(i, k - 1) = t_u + U(k, k - 1) * U(i, k);
                            }
                        }
                    }else{ /* step.31 */
                        std::swap(B_sq_norm(k), B_sq_norm(k - 1));
                        for (int i = k; i <= l - 1; i++){
                            std::swap(U(i - 1, k), (U(i - 1, k - 1)));
                        }
                        k = std::max(k - 1, 1);
                    }
                }
            }
            if(!FLAG_startagain){
                g++;
            }
        }
    }
}