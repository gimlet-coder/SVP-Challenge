#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> 
#include <stdexcept> 
#include <iomanip> //  表示の整形用
#include <chrono>  //  時間計測用

#include "lattice_types.hpp"

#include "DeepLLL.hpp"



// アルゴリズム 7: DeepLLL 基底簡約アルゴリズム
void DeepLLL(IntMatrix &B, const Real delta){
    if(delta <= 0.25 || 1 <= delta){
        throw std::out_of_range("Size_reduce: 不正なパラメータ δ です");
    }

    int n = B.rows();
    if(n < 2) return; 

    RealMatrix B_star, U; 
    Gram_Schmidt(B, B_star, U);

    RealVector B_norm = B_star.rowwise().squaredNorm(); // B_i = ||b*_i||^2 step.2

    int k_idx = 1;

    //進捗状況確認ツール
    auto last_print_time = std::chrono::steady_clock::now();
    long long iteration_count = 0;

    // step.4: while k <= n
    while(k_idx < n){

        // 時間をチェックして一定間隔でのみ出力する
        iteration_count++;
        auto now = std::chrono::steady_clock::now();
        
        // 100ミリ秒 (0.1秒) 以上経過していたら表示を更新
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print_time).count() > 100) {
            
            // 進捗バーの作成 (現在の k_idx の位置を視覚化)
            int bar_width = 20; 
            int pos = (k_idx * bar_width) / n;
            std::string progress_bar = "[" + std::string(pos, '#') + std::string(bar_width - pos, ' ') + "]";

            // 第1ベクトルの長さ(ノルム)を表示用にdoubleにキャスト
            // (多倍長型のままだとstd::scientificなどが効かない場合がある)
            double current_norm = 0.0;
            try {
                // ルートを取って長さにする
                current_norm = static_cast<double>(boost::multiprecision::sqrt(B.row(0).squaredNorm()));
            } catch (...) {
                current_norm = -1.0; // エラー時は-1表示
            }

            // std::cerr に出力 (標準出力 std::cout は計算結果用にとっておく)
            std::cerr << "\r" // 行頭に戻る
                      << "Iter: " << std::setw(8) << iteration_count 
                      << " | k: " << std::setw(3) << k_idx << "/" << n 
                      << " " << progress_bar 
                      << " | ||b_0||: " << std::scientific << std::setprecision(4) << current_norm
                      << std::flush; // バッファを強制出力して即座に表示

            last_print_time = now;
        }
        
        for (int j = k_idx - 1; j >= 0; j--){
            Size_reduce_partial(B, U, k_idx, j);// Step 6: Size-reduce(k, j)
        }
              
        Real C_scalar = B.row(k_idx).cast<Real>().squaredNorm();
        int i_idx = 0;

        while(i_idx < k_idx){       
            if(C_scalar >= delta * B_norm(i_idx)){
                C_scalar -= U(k_idx, i_idx) * U(k_idx, i_idx) * B_norm(i_idx);
                i_idx++;
            }else{
                IntVector temp = B.row(k_idx);
                for (int j = k_idx; j > i_idx; j--){
                    B.row(j) = B.row(j - 1);
                }
                B.row(i_idx) = temp;
                GSOUpdate_DeepLLL_partial(U, B_norm, i_idx + 1, k_idx + 1);
                k_idx = std::max(i_idx, 1) - 1;
                
            }
        }
        k_idx++;
    }
    std::cerr << std::endl;
}