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
    int h = B.rows(), z = h - 1, g = 1;
    Matrix B_star, U; // GSOを求める際の受け皿
    Gram_Schmidt(B, B_star, U);
    Vector B_sq_norm(h); // 疑似コードでいうところの B (B_i = ||b^*_i||^2) にあたるもの
    for (int i = 0; i < h; i++){
        B_sq_norm(i) = B_star.row(i).squaredNorm(); // GSOベクトルの2乗ノルムが B_sq_norm(i) に入るよう調整
    }
    
    while(g <= z){
        if(B_sq_norm(g) < 1e-20){ // ノルムの2乗が十分小さければ0ベクトルであると判定する(浮動小数点の誤差対策)
            // Step4-9: B_g=0 なら末尾と交換し z-- (0ベクトル除去)
            if(g < z){
                B.row(g).swap(B.row(z));
            }
            z--;
            Gram_Schmidt(B, B_star, U);
            for (int i = 0; i < h; i++){
                B_sq_norm(i) = B_star.row(i).squaredNorm(); // GSO情報を更新するたびに B も更新することを忘れないように注意
            }
        }
        int l = g, k = g;
        bool FLAG_startagain = false;
        while(k <= l && !FLAG_startagain){
            Size_reduce_partial(B, U, k, k - 1); //step.14
            Gram_Schmidt(B, B_star, U);
            for (int i = 0; i < h; i++) {
                B_sq_norm(i) = B_star.row(i).squaredNorm();
            }
            Scalar nu = U(k, k - 1), B_proj = B_sq_norm(k) + nu*nu * B_sq_norm(k - 1); //step.13
            if(B_proj >= delta * B_sq_norm(k - 1)){ // step.15 Lovasz 条件
                for(int j = k - 2; j >= 0; j--){
                    Size_reduce_partial(B, U, k, j);
                }
                // Size-reduce 後の B_k の再計算 (B_norm(k_idx) の更新)
                Gram_Schmidt(B, B_star, U);
                for (int i = 0; i < h; i++){
                    B_sq_norm(i) = B_star.row(i).squaredNorm(); // GSOベクトルの2乗ノルムが B_sq_norm(i) に入るよう調整
                }
                k++;
            }else{ //step.16
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
                    Scalar nu = U(k, k - 1);
                    Scalar B_proj = B_sq_norm(k) + nu * nu * B_sq_norm(k - 1);

                    // nu および B の更新
                    Scalar t = B_sq_norm(k - 1) / B_proj;
                    Scalar inv_nu = Scalar(1) / nu; // ν = μ_{k,k-1}
                    Scalar nu_old = nu; // 退避
                    // step. 21
                    B.row(k).swap(B.row(k - 1));

                    // B の更新
                    B_sq_norm(k) = B_proj;
                    B_sq_norm(k - 1) = t * B_sq_norm(k - 1);

                    // μ の更新
                    // まず μ_{k, j} (j <= k-2) と μ_{k-1, j} の入れ替え・更新
                    //step. 22
                    for (int j = 0; j <= k - 2; j++) {
                        Scalar mu_kj_old   = U(k, j);
                        Scalar mu_k1j_old  = U(k - 1, j);
                        // μ'_{k-1,j} = mu_kj_old
                        U(k - 1, j) = mu_kj_old;
                        // μ'_{k,j}   = mu_k1j_old - nu * mu_kj_old
                        U(k, j) = mu_k1j_old - nu_old * mu_kj_old;
                    }

                    // μ_{k, k-1} と μ_{k-1, k-1}=1 の更新
                    U(k, k - 1) = nu_old * t;   // μ'_{k, k-1} = ν * t
                    // μ_{k-1, k-1} は 1 のまま

                    // μ_{i, k-1}, μ_{i, k} (i > k) の更新
                    for (int i = k + 1; i <= l; i++) {
                        Scalar mu_ik_old  = U(i, k);
                        Scalar mu_ik1_old = U(i, k - 1);
                        // μ'_{i, k-1} = mu_ik_old + nu * mu_ik1_old
                        U(i, k - 1) = mu_ik_old + nu_old * mu_ik1_old;
                        // μ'_{i, k}   = mu_ik1_old - nu * mu_ik_old
                        U(i, k) = mu_ik1_old - nu_old * mu_ik_old;
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
