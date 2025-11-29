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
    int h = B.rows(), n = B.cols(), z = h - 1, g = 0;
    Matrix B_star(h, n), U(h, h); // GSOを求める際の受け皿
    U.setZero();
    Vector B_sq_norm(h); // 疑似コードでいうところの B (B_i = ||b^*_i||^2) にあたるもの
    bool FLAG_startagain = false;
    while(g <= z){
        //step.5
        B_star.row(g) = B.row(g);
        for (int j = 0; j < g; j++){
            if(B_sq_norm(j) > Scalar("1e-20")){// ゼロ除算対策
                U(g, j) = B.row(g).dot(B_star.row(j)) / B_sq_norm(j);
            }else{
                U(g, j) = 0;
            }
            B_star.row(g) -= U(g, j) * B_star.row(j);
        }
        B_sq_norm(g) = B_star.row(g).squaredNorm();          
        U(g, g) = 1.0;

        
        if(g == 0){
            g = 1;
            continue;
        }else{
            //step.11
            int l = g, k = g; // l, k は 0-based index
            FLAG_startagain = false;
            while(k <= l && !FLAG_startagain){
                Size_reduce_partial(B, U, k, k - 1); //step.14
                Scalar nu = U(k, k - 1);
                Scalar B_proj = B_sq_norm(k) + nu * nu * B_sq_norm(k - 1);


                if(B_proj >= delta * B_sq_norm(k - 1) - 1e-20){ // step.15 Lovasz 条件
                    for(int j = k - 2; j >= 0; j--){
                        Size_reduce_partial(B, U, k, j);
                    }
                    k++;
                }else{ //step.16
                    /* step.17 */
                    if(B.row(k).squaredNorm() < Scalar("1e-20")){ //上記同様浮動小数の誤差対策で判定基準を設ける
                        if(k < z){                            
                            B.row(k).swap(B.row(z));
                        }
                    z--;
                    g = k;
                    FLAG_startagain = true;
                    }else{
                        // step. 21
                        B.row(k).swap(B.row(k - 1));
                        //step. 22
                        for (int j = 0; j <= k - 2; j++) {
                            U(k, j).swap(U(k - 1, j));
                        }
                        if(B_proj > Scalar("1e-20")){
                            if(B_sq_norm(k) < Scalar("1e-20")){
                                B_sq_norm(k - 1) = B_proj;
                                B_star.row(k - 1) *= nu;
                                U(k, k - 1) = 1.0 / nu;
                                for (int i = k + 1; i <= l; i++){
                                    U(i, k - 1) /= nu;
                                }
                            }else{
                                //step.28
                                Scalar t = B_sq_norm(k - 1) / B_proj;
                                U(k, k - 1) = nu * t;
                                Vector w = B_star.row(k - 1);
                                B_star.row(k - 1) = B_star.row(k) + nu * w.transpose();
                                B_sq_norm(k - 1) = B_proj;
                                if(k <= l){
                                    B_star.row(k) = -U(k, k - 1) * B_star.row(k) + (B_sq_norm(k) / B_proj) * w.transpose();
                                    B_sq_norm(k) *= t;
                                }// ここまで step.28
                                for (int i = k + 1; i <= l; i++){
                                    Scalar temp = U(i, k);
                                    U(i, k) = U(i, k - 1) - nu * temp;
                                    U(i, k - 1) = temp + U(k, k - 1) * U(i, k);
                                }
                            } // step.30
                        }else{
                            B_sq_norm(k).swap(B_sq_norm(k - 1));
                            B_star.row(k).swap(B_star.row(k - 1));
                            for (int i = k + 1; i <= l; i++){
                                U(i, k).swap(U(i, k - 1));
                            }
                        }
                        // k を戻す
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