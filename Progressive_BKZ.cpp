
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <stdexcept>

#include <gmp.h>
#ifdef __linux__
#include <fplll/bkz.h>
#include <fplll/lll.h>
#endif

#include <boost/multiprecision/cpp_int.hpp>

#include "lattice_types.hpp"
#include "Progressive_BKZ.hpp"

#ifdef __linux__
/*----------------------------------------------------------------------------
 Eigen <-> fplll 変換ヘルパ
/ ----------------------------------------------------------------------------
*/
static void eigen_to_fplll(const IntMatrix &B, fplll::ZZ_mat<mpz_t> &Z) {
    int rows = B.rows();
    int cols = B.cols();
    Z.resize(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Integer (cpp_int) -> string -> mpz_t
            std::string v = B(i, j).str();
            if (mpz_set_str(Z(i, j).get_data(), v.c_str(), 10) != 0) {
                throw std::runtime_error("mpz_set_str failed");
            }
        }
    }
}

static void fplll_to_eigen(const fplll::ZZ_mat<mpz_t> &Z, IntMatrix &B) {
    int rows = Z.get_rows();
    int cols = Z.get_cols();
    B.resize(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            char *s = mpz_get_str(nullptr, 10, Z(i, j).get_data());
            B(i, j) = Integer(s);
            free(s);
        }
    }
}

static bool run_fplll_bkz(IntMatrix &B, int beta, double delta) {
    fplll::ZZ_mat<mpz_t> Z;
    eigen_to_fplll(B, Z);

    std::vector<fplll::Strategy> strategies;
    try {
        std::string strat_path = fplll::strategy_full_path(fplll::default_strategy());
        strategies = fplll::load_strategies_json(strat_path);
    } catch (...) {
        // 戦略ファイルが見つからない場合は空の戦略を用意
        strategies.emplace_back(fplll::Strategy::EmptyStrategy(beta));
    }

    fplll::BKZParam params(beta, strategies, delta, fplll::BKZ_DEFAULT, 1);
    int rc = fplll::bkz_reduction(&Z, nullptr, params);
    if (rc != 0) {
        return false;
    }

    fplll_to_eigen(Z, B);
    return true;
}
#endif

/* ----------------------------------------------------------------------------
 基底のランダム化 (Unimodular変換) //あえて行列をほどよく汚して再度簡約させてみる
/ ----------------------------------------------------------------------------
*/
static void randomize_basis(IntMatrix &B, int num_operations) {
    int n = B.rows();
    std::random_device rd;
    std::mt19937 gen(rd()); // 適当な乱数を生成
    std::uniform_int_distribution<> dist_idx(0, n - 1);
    std::uniform_int_distribution<> dist_sign(0, 1);

    for (int k = 0; k < num_operations; k++) {
        int i = dist_idx(gen);
        int j = dist_idx(gen);
        
        if (i == j) continue;

        // 行操作: row(i) = row(i) ± row(j)
        if (dist_sign(gen)) {
            B.row(i) += B.row(j);
        } else {
            B.row(i) -= B.row(j);
        }
    }
}
/*----------------------------------------------------------------------------
ファイルを出力する関数
/ ----------------------------------------------------------------------------
*/


static void save_current_state(const IntMatrix &B, int n, int beta, const Real best_norm_sq, const std::string &prefix) {
    // ファイル名: [prefix]_beta_[beta].txt の形式
    // Global Best Record の場合は beta の代わりに "best" を使う
    std::string beta_str = (beta == -1) ? "best" : std::to_string(beta);
    std::string filename = prefix + "_beta_" + beta_str + ".txt";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "警告: 状態ファイル " << filename << " を開けませんでした。" << std::endl;
        return;
    }

    // ノルムは二乗ノルムとして格納 (多倍長 Real 型)
    Real b1_norm_sq = B.row(0).cast<Real>().squaredNorm();
    double current_norm = std::sqrt(static_cast<double>(b1_norm_sq));

    outfile << "--- Progressive BKZ State Save ---" << std::endl;
    outfile << "Dimension: " << n << std::endl;
    outfile << "Last Completed Beta: " << beta_str << std::endl;
    outfile << "Shortest Norm Squared (||b1||^2): " << std::setprecision(30) << b1_norm_sq << std::endl;
    outfile << "Shortest Vector Norm (||b1||): " << std::fixed << std::setprecision(4) << current_norm << std::endl;
    outfile << "Best Norm Found in This Run: " << std::fixed << std::setprecision(4) << std::sqrt(static_cast<double>(best_norm_sq)) << std::endl;
    outfile << "--- Basis Matrix B (n x n) ---" << std::endl;
    
    // 基底行列 B の出力
    for (int i = 0; i < n; i++) {
        // 各行ベクトルを [] で囲んで出力（読み込み関数が処理しやすいように）
        outfile << "[";
        for (int j = 0; j < n; j++) {
            // 多倍長浮動小数点から整数への丸め込み
            long long component = static_cast<long long>(B(i, j));
            outfile << component << (j < n - 1 ? " " : "");
        }
        outfile << "]" << std::endl;
    }
    
    outfile.close();
    std::cout << "\n[SAVE] 状態を " << filename << " に保存しました。ノルム: " 
              << std::fixed << std::setprecision(4) << current_norm << std::endl;
}

/*----------------------------------------------------------------------------
 データの読み込み
/ ----------------------------------------------------------------------------
*/
static void load_data(const std::string &filename, int n, IntMatrix &B){
    if (n <= 0) {
        throw std::invalid_argument("次元 n は 1 以上である必要があります。");
    }
    // scanning_lattice.cpp の機能を使って、ファイルから全てのデータを読み込む
    // n x n の行列としてファイルから読み込みをする
    B = load_challenge_matrix(filename, n, n);

    
    std::cout << "ファイル読み込みが完了しました。次元: " << n << std::endl;
}

/* ----------------------------------------------------------------------------
 Randomized Progressive BKZ 
@brief 最短ベクトルが求まったとしても, 基底をある程度汚くして max_retries 回実行して最短であるかを確認する
 ----------------------------------------------------------------------------
*/

// グローバル保持することによって途中で終了しても再開できるようにする
static Real Global_best_norm_sq("1e100"); // 十分大きな値
static IntMatrix Global_B_best;


void Progressive_BKZ(const std::string &filename, int n, int start_beta, int max_beta, Real delta, double target_norm, int max_retries) {
    IntMatrix B(n, n);
    
    // 1. データ読み込み
    load_data(filename, n, B);
    
#if defined(__linux__)
    // 1.5 fplll で前処理 BKZ（整数BKZで軽く整形）
    int preprocess_beta = std::min(n, std::max(start_beta, 20));
    double preprocess_delta = static_cast<double>(delta);
    std::cout << "[fplll] BKZ preprocess (beta=" << preprocess_beta << ", delta=" << preprocess_delta << ")..." << std::endl;
    if (run_fplll_bkz(B, preprocess_beta, preprocess_delta)) {
        std::cout << "[fplll] preprocess done. ||b1|| = "
                  << std::sqrt(static_cast<double>(B.row(0).squaredNorm())) << std::endl;
    } else {
        std::cerr << "[fplll] preprocess failed. 続行しますが前処理なしになります。" << std::endl;
    }
#endif



    Global_best_norm_sq = B.row(0).cast<Real>().squaredNorm();
    Global_B_best = B;

    std::cout << "--- Initial State ---" << std::endl;
    std::cout << "||b_1||^2 = " << Global_best_norm_sq << std::endl;

    std::cout << "\n=== Start Progressive BKZ ===" << std::endl;
    // 2. パラメータ設定

    auto total_start = std::chrono::high_resolution_clock::now();

    std::cout << "\n=== Start Progressive BKZ (Beta: " << start_beta << " -> " << max_beta << ") ===" << std::endl;

    for (int beta = start_beta; beta <= max_beta; ) { // POINT!! ベータはループの最下層で増加させる


        std::cout << "\n----------------------------------------" << std::endl;
        std::cout << "[Phase Beta = " << beta << "]" << std::endl;
        
        Real best_norm_sq = B.row(0).cast<Real>().squaredNorm();
        int retry_count = 0;
        bool improved_in_phase = false;

        while (retry_count < max_retries) {
            auto iter_start = std::chrono::high_resolution_clock::now();

            // BKZ関数を呼び出し
            DeepBKZ(B, beta, delta);

            auto iter_end = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end - iter_start).count() / 1000.0;

            Real current_sq = B.row(0).cast<Real>().squaredNorm();
            double current_norm = std::sqrt(static_cast<double>(current_sq));

            // 進捗表示
            std::cout << "  Iter " << std::setw(2) << retry_count + 1 
                      << " | Norm: " << std::fixed << std::setprecision(4) << current_norm 
                      << " | Time: " << std::setprecision(2) << elapsed << "s";

            // 目標達成チェック
            if (current_norm <= target_norm) {
                std::cout << "\n\n[SUCCESS] 目標到達！ (Norm: " << current_norm << ")" << std::endl;
                // ここで、初期設定していた目標値に達したことを表している。計算は続行する
            }

            // 改善判定
            if (current_sq < best_norm_sq - 1e-5) { // ここの精度はざっくりでいいので10^{-5} 程度にした

                best_norm_sq = current_sq;
                Global_best_norm_sq = current_sq;
                Global_B_best = B;
                save_current_state(Global_B_best, n, -1, Global_best_norm_sq, "lattice_state_80d"); // ここの 80 は次元を変えるたびに変更したほうがいい
            std::cout << "\n!!! NEW GLOBAL RECORD FOUND: " << std::sqrt(static_cast<double>(Global_best_norm_sq)) << " !!!\n";
            
                best_norm_sq = current_sq;
                retry_count = 0; // 改善したらカウントリセット
                improved_in_phase = true;
            } else {
                std::cout << " -> Stalled." << std::endl;
                retry_count++;
                
                // 停滞時: ランダム化して再挑戦
                if (retry_count < max_retries) {
                    std::cout << "    -> Randomizing basis..." << std::endl;
                    // ベクトルをかき混ぜる (次元数に応じて回数を調整)
                    randomize_basis(B, n * 2); 
                }
            }
        }

        if (beta < max_beta) {
        int next_beta = beta + 1; // 少なくとも + 1

        if (beta >= 30) {
            next_beta = beta + 5;
        } else if (beta >= 25) {
            next_beta = beta + 3;
        }
        
        // 最終チェック: max_beta を超えないようにする
        if (next_beta > max_beta) {
            beta = max_beta; // 次の周は max_beta
        } else {
            beta = next_beta; // 通常の増分
        }
        // ここでは for の制御部で beta++ しないため、beta = next_beta で次の周の beta が決定する
        
        } else {
            // beta == max_beta の実行が終わったら、ループを終了させる
            break;
        }

        B = Global_B_best;

    }


    auto total_end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start).count();

    // 3. 最終結果
    std::cout << "\n========================================" << std::endl;
    std::cout << "全工程完了 (Total Time: " << total_time << " s)" << std::endl;
    std::cout << "最終最短ベクトルノルム: " << std::sqrt(static_cast<double>(B.row(0).squaredNorm())) << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "最終基底 B (Top 5 rows):\n" << B.topRows(5) << std::endl; // おまけで短いベクトルを5つ列挙する
}
