#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <iomanip> //  表示の整形用
#include <chrono>  //  時間計測用


#include "lattice_types.hpp"

bool ENUM(const Matrix &U, const Vector &B_norm, const Vector &R_squares, Vector &v_out, const int k_begin, const int k_end, long long &node_count){
    const int n_rows = static_cast<int>(U.rows());
    const int n_cols = static_cast<int>(U.cols());

    if (k_begin > k_end) throw std::out_of_range("k_begin と k_end の大小関係が不正です");
    if (k_begin < 0) throw std::out_of_range("k_begin は 0 以上である必要があります");
    if (n_rows <= k_end) throw std::out_of_range("k_end が行数の範囲外です");
    if (n_rows != n_cols || n_rows != static_cast<int>(B_norm.rows())) {
        throw std::out_of_range("GSO 行列 U と ノルム B_norm のサイズが一致しません");
    }

    /* 初期チェックと準備 */

    const int n = k_end - k_begin + 1;

    // Sigma, rho, r のサイズ準備 (インデックスエラー防止のため少し余裕を持つ)
    Matrix Sigma(n + 2, n + 2);
    Sigma.setZero();
    Vector rho(n + 2);
    rho.setZero();
    std::vector<int> r(n + 2); // step.2,3: r_0 = 0, r_i = i
    r[0] = 0;
    for (int i = 1; i <= n; i++){
        r[i] = i;
    }
    bool FLAG_found = false;
    node_count = 0;

    // v_temp, c, w_int の初期化
    Eigen::Matrix<long long, Eigen::Dynamic, 1> v_temp(n);
    v_temp.setZero();
    Vector c(n);
    c.setZero();
    Eigen::Matrix<long long, Eigen::Dynamic, 1> w_int(n);
    w_int.setZero();

    int last_nonzero = 1; // step.7
    int k = n;            // step.8: 根 (k を増やして葉へ)

    // --- 進捗表示用の変数 ---
    auto start_time = std::chrono::steady_clock::now();
    auto last_print_time = start_time;

    while(true){ // step.9
        node_count++;
        // --- 進捗表示ロジック ---
        // 2^16回に1回チェックする
        if ((node_count & 0xFFFF) == 0) { 
            auto now = std::chrono::steady_clock::now();
            // 前回の表示から0.5秒以上経過していたら出力
            if (std::chrono::duration_cast<std::chrono::milliseconds>(now - last_print_time).count() > 500) {
                double elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
                double nps = node_count / (elapsed + 1e-9); // nodes per second

                std::cerr << "\r" // 同じ行を上書き
                          << "[ENUM] Nodes: " << std::setw(10) << node_count 
                          << " | Depth(k): " << std::setw(2) << k // 現在の探索深さ
                          << " | Speed: " << std::scientific << std::setprecision(2) << nps << " n/s"
                          << std::flush;
                
                last_print_time = now;
            }
        }


        const int global_k = k + k_begin - 1;
        Scalar diff = static_cast<Scalar>(v_temp(k - 1)) - c(k - 1);

        rho(k) = rho(k + 1) + diff * diff * B_norm(global_k); // step.10


        if(rho(k) <= R_squares(k - 1)){ // step.11
            if (k == 1) { // 葉 k = 1
                // 0 ベクトルは捨てて次候補へ
                bool all_zero = true;
                for (int i = 0; i < n; i++) {
                    if (v_temp(i) != 0) {
                        all_zero = false;
                        break;
                    }
                }
                if (!all_zero) { // 非零ベクトル発見
                    v_out = v_temp.cast<Scalar>();
std::cerr << "\n見つかったノルム: " << R_squares(k - 1) << " (Nodes: " << node_count << ")" << std::endl;
                    return true;
                }
                // 0 ベクトルなら次の候補へ step.33 ~ 38 と同じ処理を行う
                if (v_temp(k - 1) > c(k - 1)){
                    v_temp(k - 1) += w_int(k - 1);
                }else{
                    v_temp(k - 1) -= w_int(k - 1);
                }
                w_int(k - 1)++;
                continue;
            }
            k--; // step.15:
            r[k] = std::max(r[k], r[k + 1]); // step.16

            // step.17-19: Sigma の更新（添字を0-basedに揃える）
            for (int i_local = r[k]; i_local >= k + 1; i_local--) {
                int g_i = i_local + k_begin - 1;
                int g_k = k + k_begin - 1;
                Sigma(i_local, k - 1) = Sigma(i_local + 1, k - 1) + static_cast<Scalar>(v_temp(i_local - 1)) * U(g_i, g_k);
            }
            // step.20-23: c_k, v_k, w_k の設定
            c(k - 1) = -Sigma(k + 1, k - 1);
            Scalar rounded_c_k = boost::multiprecision::round(c(k - 1));
            v_temp(k - 1) = static_cast<long long>(rounded_c_k);
            w_int(k - 1) = 1;
        } else { // rho(k) > R_square
            k++; // step.24
            if (k == n + 1) { // step.25
                return FLAG_found; // step.26: 目的の v \in L が存在しない場合の終了
            }

            r[k - 1] = k; // step.28

            if (k >= last_nonzero){ // step.29
                last_nonzero = k;
                v_temp(k - 1) += 1; // step.31: v_k <- v_k + 1
                w_int(k - 1) = 1; // w をリセット
            }else{
                if (v_temp(k - 1) > c(k - 1)){ // step.32-38: 交互にずらしながら探索
                v_temp(k - 1) -= w_int(k - 1); // step.34
                }else{
                    v_temp(k - 1) += w_int(k - 1); // step.36
                }
            w_int(k - 1)++; // step.37
            }
        }
    }
}
