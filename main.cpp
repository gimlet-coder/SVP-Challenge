#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <cmath> // std::fabs() (絶対値の計算) std::round() に用いる 
#include <stdexcept> // std::out_of_range のために必要

#include <random> // 10^7 以上の疑似乱数を生成するために必要
#ifdef _WIN32
#include <windows.h> // 文字化け対策 (Windows 専用)
#endif
#include <iomanip> // 表示桁数を確保するのに使う std::fixed << std::setprecision(0)
#include <chrono> // 計算時間計測するために使う

#include "lattice_types.hpp"
// 上記のヘッダーによって Vector と Matrix は Scalar 型になっている



#include "main.hpp"

const std::string FILE_NAME = "lattice.txt"; // スクリプトが出力するファイル名と合わせる

const int DIM = 40; // 次元を入力

const int MAX_BETA = 40;
const Real DELTA = 0.999;
const double TARGET_NORM = 2100.0; // 64次元の一番下の記録よりは超えたい

#if 1

int main(int argc, char* argv[]){
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    int start_beta = 40; // コマンドライン引数がない場合のデフォルト値
    int max_retries = 5;

    if (argc > 1){
        start_beta = std::stoi(argv[1]);
    }
    if (argc > 2){
        max_retries = std::stoi(argv[2]);
    }

    try {
        Progressive_BKZ(FILE_NAME, DIM, start_beta, MAX_BETA, DELTA, TARGET_NORM, max_retries);
    } catch (const std::exception& e) {
        std::cerr << "エラーが発生しました: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

#endif


#if 0
// BKZの動作確認
int main(){
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    int n = 40; // 次元を指定
    Matrix B(n, n);
    // lattice_file.txt から行列を読み込みたい
    std::string filename = "lattice_file.txt";
    std::cout << "以下のファイルから行列を取得中 --- " << filename << " ---" << std::endl;
    try{

        Matrix B = load_challenge_matrix(filename, n, n);
        std::cout << "--- Initial Basis (First row) ---" << std::endl;
        // うまく読み込みできているか最初の行で呼び出しをかけてみる
        std::cout << B.row(0) << std::endl; 
        std::cout << "||b_1||^2 = " << B.row(0).squaredNorm() << std::endl;

        // BKZのパラメータ設定
        // n = 40 なので beta = 40 が最強設定
        int beta = 40; 
        Scalar delta = 0.99;

        std::cout << "\nBKZ(beta=" << beta << ", delta=" << delta << ") を実行中..." << std::endl;

        auto start_time = std::chrono::high_resolution_clock::now(); // ここから計測を始める
        BKZ(B, beta, delta);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        std::cout << "\n--- BKZ実行後 ---" << std::endl;
        Scalar b1_sq = B.row(0).squaredNorm();
        std::cout << "||b_1||^2 = " << b1_sq << std::endl;
        std::cout << "\n--- BKZ 実行時間 ---\n";
        std::cout << "実行時間: " << duration.count() << " ms (" << static_cast<double>(duration.count()) / 1000.0 << " s)\n";
        std::cout << "------------------------\n\n";
        std::cout << "\n簡約基底 B:\n" << B << std::endl;
    }catch (const std::exception& e) {
        std::cerr << "エラーが発生しました: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

#endif


#if 0

int main() {
#ifdef _WIN32
    SetConsoleOutputCP(65001);
#endif
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

    // --- 例 3.2.1 の基底行列 B ---
    int n = 10;
    Matrix B(n, n);
    B <<  -79, 35,  31,  83,  -66, 35, -32, 46, 21, 2,
          43, -64, -37, -31, -27, -7, -42, 21, 16, 16,
          -1, -97, -91, -43, 19, -21, -65, -36, 34, -55,
          -58, -38, 87, 42, 94, -83, 66, -69, -2, -30,
          84, -61, 93, -67, 3, 94, 31, 27, -60, 98,
          -1, 34, 58, -38, 29, 67, -18, 15, -75, -16,
          19, 16, 52, 32, -20, 55, 94, -34, 4, 80,
          -58, -17, 99, 93, -49, -53, 24, 51, 5, 93,
          17, 31, 78, 53, 40, -22, -39, 7, 70, -98,
          93, -6, -7, -12, 79, -40, 27, -95, 98, 20;

    std::cout << "入力した基底行列B:\n" << B << std::endl << std::endl;
/*
    Scalar delta = 0.26;
    std::cout <<"delta = " << delta << " におけるDeepLLLの実行中" << std::endl;
    DeepLLL(B, delta);
    std::cout << "DeepLLLの実行結果:\n" << B << "\n\n";
*/
    // --- GSOの計算 ---
    Matrix B_star(n, n);
    Matrix U(n, n);
    Gram_Schmidt(B, B_star, U);

    // ENUM用にノルム配列を作成
    Vector B_norm(n);
    for(int i = 0; i < n; i++){
        B_norm(i) = B_star.row(i).squaredNorm();
    }
    std::cout << "GSOの2乗ノルム (B*):\n" << B_norm.transpose() << std::endl << std::endl;

    // --- 探索設定 ---
    // 第1基底ベクトル b1 のノルム
    Scalar b1_sq_norm = B.row(0).squaredNorm();
    std::cout << "||b1||^2 = " << b1_sq_norm << "\n";

    // b1 より短いベクトルを探す (R^2 = ||b1||^2 - 1)
    Scalar R_square = b1_sq_norm - 1; 
    
    std::cout << "||v||^2 <= " << R_square << " となるような v を探索中 ...\n";

    Vector v_coeffs(n);
    long long node_count = 0;
    bool found = ENUM(U, B_norm, R_square, v_coeffs, 0, n - 1, node_count);

    if (found) {
        std::cout << "\n!!! 最短ベクトルを発見 !!!\n";
        std::cout << "係数 (v1...v10): " << v_coeffs.transpose() << "\n";

        // 実際の格子ベクトルを計算
        Vector v_lattice = Vector::Zero(n);
        for(int i = 0; i < n; i++){
            // v_coeffs(i) * B.row(i)
            v_lattice += v_coeffs(i) * B.row(i);
        }
        std::cout << "格子ベクトル: " << v_lattice.transpose() << "\n";
        std::cout << "2乗ノルム: " << v_lattice.squaredNorm() << "\n";

        // 教科書の正解チェック (0, 1, 1, 0, 1) or (0, -1, -1, 0, -1)
    } else {
        std::cout << "\n最短ベクトルは見つかりませんでした。b1 が最短だと推測されます。\n b1 = " << B.row(0) <<std::endl;
    }

    return 0;
}
#endif
