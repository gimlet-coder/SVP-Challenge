#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include "BKZ.hpp"
#include "lattice_types.hpp" 

// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている
// MULTI_PRECISION_EPSILON は上記ヘッダーにより定義されている この値を用いて非ゼロ判定をする

/*memo:
B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す
hoge.dot(huga) は hoge と huga の内積を表す
hoge.squaredNorm() はベクトルの長さの2乗を表す
*/

void BKZ (Matrix &B, int beta, const Scalar delta){
    if(delta < 0.25 || 1 <= delta){
        throw std::out_of_range("Size-reduce 不正なインデックス delta です (1/4 < delta < 1)");
    }
    if(beta < 2 || beta > B.rows()){
        throw std::out_of_range("Block size 不正なインデックス beta です (2 <= beta <= n)");
    }
    /* --- 引数チェック完了 */
    int n = B.rows(); // 次元を n とする
    LLL(B, delta); // step.2 LLL簡約で基底を更新する

    //ここで, LLL簡約後の GSO 情報を保存したい

    Matrix B_star, U; // GSO 情報の受け皿
    Gram_Schmidt(B, B_star, U);

    Vector B_norm(n);
    for(int i = 0; i < n; i++){
        B_norm(i) = B_star.row(i).squaredNorm();
    }

    Scalar b1_sq_norm = B.row(0).squaredNorm();
    int z = 0, k = -1; // ここで, k を 0-based に変更する z はあくまでカウンターなので 0 でOK

    while(z < n - 1){
        k = (k + 1) % (n - 1);
        int l = std::min(k + beta - 1, n - 1);
        int h = std::min(l + 1, n);
        Scalar R_square = b1_sq_norm * 0.99;
        Vector v_coeffs; // 部分射影格子 L 上の最短な非ゼロベクトルの係数ベクトルの保存先
        long long node_count = 0; // ENUMでどのくらいのノード数を要したか保存する
        bool found = ENUM(U, B_norm ,R_square ,v_coeffs , k, l, node_count);

        
        if(found){
            Vector v_lattice = Vector::Zero(B.cols());
            for (int i = 0; i < v_coeffs.size(); i++){
                v_lattice += static_cast<Scalar>(v_coeffs(i)) * B.row(k + i);
            }
            Scalar proj_v_sq_norm = v_lattice.squaredNorm();
            for (int i = 0; i < k; i++){
                Scalar mu_v_j = v_lattice.dot(B_star.row(i)) / B_norm(i); // <v, b*_j> / ||B_i||
                proj_v_sq_norm -= mu_v_j * mu_v_j * B_norm(i); // 射影成分を除く
            }
            if(proj_v_sq_norm < MULTI_PRECISION_EPSILON){ // 計算誤差分よりも小さい場合は 0 とする
                proj_v_sq_norm = 0;
            }
            
            if(B_norm(k) - proj_v_sq_norm > MULTI_PRECISION_EPSILON){
                z = 0; // カウンターを 0 にセット

                // MLLL を呼び出すための準備
                Matrix B_h_added_v(h + 1, B.cols()); // b_0 ~ b_{h-1} に加えて v があるので h + 1 個の要素がある
                int idx = 0;
                for (int i = 0; i < h + 1; i++){
                    if(i == k){
                        B_h_added_v.row(i) = v_lattice; // i = k となるところでは, idx を動かさない
                    }else{
                        B_h_added_v.row(i) = B.row(idx);
                        idx++;
                    }
                }            
                MLLL(B_h_added_v, delta); // この段階で, B_h_added_v は (h + 1) 行あるが, 最後の行は 0 ベクトルなので消す
                for (int i = 0; i < h; i++){ // B_h_added_v.row(h) は 0 ベクトル
                    B.row(i) = B_h_added_v.row(i);
                }
                Gram_Schmidt(B, B_star, U); // OPTIMIZE: 最終的にはGSOUpdate-BKZを実装して差し替える
            }else{
                goto NO_UPDATE;
            }
        }else{
            NO_UPDATE:
            z++;
            LLL(B, delta);
            Gram_Schmidt(B, B_star, U);
            for(int i = 0; i < n; i++){
                B_norm(i) = B_star.row(i).squaredNorm();
            }
        }
    }
}