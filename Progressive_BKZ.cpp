
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <stdexcept>

#include "lattice_types.hpp"
#include "Progressive_BKZ.hpp"

/* ----------------------------------------------------------------------------
 基底のランダム化 (Unimodular変換) //あえて行列をほどよく汚して再度簡約させてみる
/ ----------------------------------------------------------------------------
*/
static void randomize_basis(Matrix &B, int num_operations) {
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


static void save_current_state(const Matrix &B, int n, int beta, const Scalar best_norm_sq, const std::string &prefix) {
    // ファイル名: [prefix]_beta_[beta].txt の形式
    // Global Best Record の場合は beta の代わりに "best" を使う
    std::string beta_str = (beta == -1) ? "best" : std::to_string(beta);
    std::string filename = prefix + "_beta_" + beta_str + ".txt";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "警告: 状態ファイル " << filename << " を開けませんでした。" << std::endl;
        return;
    }

    // ノルムは二乗ノルムとして格納 (多倍長 Scalar 型)
    Scalar b1_norm_sq = B.row(0).squaredNorm();
    double current_norm = std::sqrt(static_cast<double>(b1_norm_sq));

    outfile << "--- Progressive BKZ State Save ---" << std::endl;
    outfile << "Dimension: " << n << std::endl;
    outfile << "Last Completed Beta: " << beta_str << std::endl;
    outfile << "Shortest Norm Squared (||b1||^2): " << std::setprecision(30) << b1_norm_sq << std::endl;
    outfile << "Shortest Vector Norm (||b1||): " << std::fixed << std::setprecision(4) << current_norm << std::endl;
    outfile << "Best Norm Found in This Run: " << std::fixed << std::setprecision(4) << std::sqrt(static_cast<double>(best_norm_sq)) << std::endl;
    outfile << "--- Basis Matrix B (n x n) ---" << std::endl;
    
    // 基底行列 B の出力
    for (int i = 0; i < n; ++i) {
        // 各行ベクトルを [] で囲んで出力（読み込み関数が処理しやすいように）
        outfile << "[";
        for (int j = 0; j < n; ++j) {
            // 多倍長精度を活かして正確な値をスペース区切りで出力
            outfile << B(i, j) << (j < n - 1 ? " " : "");
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
static void load_and_clean_data(const std::string &filename, int n, Matrix &B){
    if (n <= 0) {
        throw std::invalid_argument("次元 n は 1 以上である必要があります。");
    }
    // scanning_lattice.cpp の機能を使って、ファイルから全てのデータを読み込む
    // n x n の行列としてファイルから読み込みをする
    Matrix B_full = load_challenge_matrix(filename, n, n);

    // 2. ゼロベクトルと重複ベクトルのチェックと削除
    std::vector<Vector> unique_vectors;
    // Set を使って、ベクトルの重複をチェックする
    std::set<std::string> seen_vectors; 

    for (int i = 0; i < n; i++) {
        Vector current_row = B_full.row(i);
        
        // a. ゼロベクトルチェック
        if (current_row.squaredNorm() < MULTI_PRECISION_EPSILON) {
            continue; // ゼロベクトルはスキップ
        }

        // b. 重複チェック: ベクトルを行列要素を文字列に変換してSetで管理
        // 正確な比較のため、高精度なScalar型を文字列に変換して比較
        std::stringstream ss;
        ss << std::fixed << std::setprecision(20); // 高精度で出力
        for (int j = 0; j < n; j++) {
            ss << current_row(j) << ",";
        }
        std::string vector_str = ss.str();

        if (seen_vectors.find(vector_str) == seen_vectors.end()) {
            // 重複がない場合
            seen_vectors.insert(vector_str);
            unique_vectors.push_back(current_row);
        }
    }

    // 3. 必要な次元数が揃っているか確認
    if (unique_vectors.size() < static_cast<size_t>(n)) {
        std::cerr << "警告: " << filename << " から " << n << " 個の独立な基底ベクトルを読み込めませんでした。"
                  << "読み込めたのは " << unique_vectors.size() << " 個です。" << std::endl;
        // 不足分がある場合、エラーとする
        throw std::runtime_error("次元不足: 完全な基底行列を作成できませんでした。");
    }

    // 4. クリーンアップされた基底を B に格納
    B = Matrix::Zero(n, n); // 初期化
    for (int i = 0; i < n; i++) {
        B.row(i) = unique_vectors[i];
    }
    
    std::cout << "ファイル読み込みとクリーンアップが完了しました。次元: " << n << std::endl;
}

/* ----------------------------------------------------------------------------
 Randomized Progressive BKZ 
@brief 最短ベクトルが求まったとしても, 基底をある程度汚くして max_retries 回実行して最短であるかを確認する
 ----------------------------------------------------------------------------
*/

// グローバル保持することによって途中で終了しても再開できるようにする
static Scalar Global_best_norm_sq = std::numeric_limits<Scalar>::max(); // 最大値で初期化
static Matrix Global_B_best;


void Progressive_BKZ(const std::string &filename, int n, int start_beta, int max_beta, Scalar delta, double target_norm, int max_retries) {
    Matrix B(n, n);
    
    // 1. データ読み込み
    load_and_clean_data(filename, n, B);

    // 初期状態表示
    Scalar current_norm_sq = B.row(0).squaredNorm();
    std::cout << "--- Initial State ---" << std::endl;
    std::cout << "||b_1||^2 = " << current_norm_sq << std::endl;
    std::cout << "||b_1||   = " << std::sqrt(static_cast<double>(current_norm_sq)) << std::endl;

    // 2. パラメータ設定

    auto total_start = std::chrono::high_resolution_clock::now();

    std::cout << "\n=== Start Progressive BKZ (Beta: " << start_beta << " -> " << max_beta << ") ===" << std::endl;

    for (int beta = start_beta; beta <= max_beta; ) { // POINT!! ベータはループの最下層で増加させる


        std::cout << "\n----------------------------------------" << std::endl;
        std::cout << "[Phase Beta = " << beta << "]" << std::endl;
        
        Scalar best_norm_sq = B.row(0).squaredNorm();
        int retry_count = 0;
        bool improved_in_phase = false;

        while (retry_count < max_retries) {
            auto iter_start = std::chrono::high_resolution_clock::now();

            // BKZ関数を呼び出し
            DeepBKZ(B, beta, delta);

            auto iter_end = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(iter_end - iter_start).count() / 1000.0;

            Scalar current_sq = B.row(0).squaredNorm();
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
                Global_B_best = B;
                save_current_state(Global_B_best, n, -1, Global_best_norm_sq, "lattice_state_50d");
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