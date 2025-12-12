#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要
#include <iomanip> //  表示の整形用
#include <chrono>  //  時間計測用

#include "DeepBKZ.hpp"
#include "lattice_types.hpp" 

// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている
// MULTI_PRECISION_EPSILON は上記ヘッダーにより定義されている この値を用いて非ゼロ判定をする

/*memo:
B.row(i) 行列 B の i 番目の行をひとつの行ベクトルとして取り出す
B.col(j) 行列 B の j 番目の列をひとつの列ベクトルとして取り出す
hoge.dot(huga) は hoge と huga の内積を表す
hoge.squaredNorm() はベクトルの長さの2乗を表す
*/

void DeepBKZ (IntMatrix &B, int beta, const Real delta){
    if(delta < 0.25 || 1 <= delta){
        throw std::out_of_range("Size-reduce 不正なインデックス delta です (1/4 < delta < 1)");
    }
    if(beta < 2 || beta > B.rows()){
        throw std::out_of_range("Block size 不正なインデックス beta です (2 <= beta <= n)");
    }
    /* --- 引数チェック完了 */
    int n = B.rows(); // 次元を n とする
#ifdef _WIN32 // windows環境ではfplllができないのでここで軽く簡約させる
    LLL(B, 0.75); // まずは大まかにLLLで簡約する ここはざっくりでいいのでわざと delta ではなく 0.75 と少し余裕をもたせる
    DeepLLL(B, delta); // DeepLLL簡約で基底を更新する
#endif
    //ここで, LLL簡約後の GSO 情報を保存したい

    RealMatrix B_star, U; // GSO 情報の受け皿
    Gram_Schmidt(B, B_star, U);

    RealVector B_norm(n);
    for(int i = 0; i < n; i++){
        B_norm(i) = B_star.row(i).squaredNorm();
    }
    int z = 0, k = -1; // ここで, k を 0-based に変更する z はあくまでカウンターなので 0 でOK

    //進捗状況確認ツール
    auto last_print_time = std::chrono::steady_clock::now();
    long long iteration_count = 0;

    while(z < n - 1){

        // 時間をチェックして一定間隔でのみ出力する
        iteration_count++;
        auto now = std::chrono::steady_clock::now();

        // 100ミリ秒 (0.1秒) 以上経過していたら表示を更新
        if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print_time).count() > 100) {
            
            // 進捗バーの作成 (現在の k_idx の位置を視覚化)
            int bar_width = 20; 
            int pos = (z * bar_width) / n;
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
                      << " | k: " << std::setw(3) << k + 1 << "/" << n 
                      << " " << progress_bar 
                      << " | ||b_0||: " << std::scientific << std::setprecision(4) << current_norm
                      << std::flush; // バッファを強制出力して即座に表示

            last_print_time = now;
        }

        k = (k + 1) % (n - 1);
        int l = std::min(k + beta - 1, n - 1);
        int h = std::min(l + 1, n);
        Real searching_radius = B_norm(k) * 0.99; // 基準の探索半径
        int block_dim = l - k + 1; // ブロックの次元数
        RealVector pruning_bounds(block_dim); // 剪定(= 枝刈り)の境界

        for (int i = 0; i < block_dim; i++){
            double ratio = static_cast<double> (i + 1) / block_dim;
            pruning_bounds(i) = searching_radius * static_cast<Real>(ratio);
        }

        IntVector v_coeffs; // 部分射影格子 L 上の最短な非ゼロベクトルの係数ベクトルの保存先
        long long node_count = 0; // ENUMでどのくらいのノード数を要したか保存する
        bool found = ENUM_fast(U, B_norm ,pruning_bounds ,v_coeffs , k, l, node_count);

        
        if(found){
            IntVector v_lattice = IntVector::Zero(B.cols());
            for (int i = 0; i < v_coeffs.size(); i++){
                v_lattice += v_coeffs(i) * B.row(k + i);
            }
            RealVector v_lat_real = v_lattice.cast<Real>();
            Real proj_v_sq_norm = v_lat_real.squaredNorm();
            for (int i = 0; i < k; i++){
                Real mu_v_j = v_lat_real.dot(B_star.row(i)) / B_norm(i); // <v, b*_j> / ||B_i||
                proj_v_sq_norm -= mu_v_j * mu_v_j * B_norm(i); // 射影成分を除く
            }
            if(proj_v_sq_norm < MULTI_PRECISION_EPSILON){ // 計算誤差分よりも小さい場合は 0 とする
                proj_v_sq_norm = 0;
            }
            
            if(B_norm(k) - proj_v_sq_norm > MULTI_PRECISION_EPSILON){
                z = 0; // カウンターを 0 にセット

#if 0 // MLLLを回避するやり方 精度を上げてMLLLが不安定となった場合切り替える

                IntVector temp_row = B.row(n - 1);
                for (int i = n - 1; i > k; i--){
                    B.row(i) = B.row(i - 1); // 後ろにシフトさせる
                }

                B.row(k) = v_lattice; // k 番目の位置に新しい最短ベクトル v を挿入する

                DeepLLL(B, delta); // DeepLLL で先頭へ最短ベクトルを配置させる
                

#endif



#if 1 //MLLLを使う場合                
                // MLLL を呼び出すための準備
                IntMatrix B_h_added_v(h + 1, B.cols()); // b_0 ~ b_{h-1} に加えて v があるので h + 1 個の要素がある
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
#endif
                Gram_Schmidt(B, B_star, U); // OPTIMIZE: 最終的にはGSOUpdate-BKZを実装して差し替える
            }else{
                goto NO_UPDATE;
            }
        }else{
            NO_UPDATE:
            z++;
            DeepLLL(B, delta);
            Gram_Schmidt(B, B_star, U);
            for(int i = 0; i < n; i++){
                B_norm(i) = B_star.row(i).squaredNorm();
            }
        }
    }
    std::cerr <<std::endl;
}