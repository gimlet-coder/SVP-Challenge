#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <Eigen/Dense>

#include "lattice_types.hpp"

// ファイルから行列を読み込む関数
Matrix load_challenge_matrix(const std::string &filename, int rows, int cols) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("ファイルを開けませんでした: " + filename);
    }

    Matrix B(rows, cols);
    std::string line;
    int row_idx = 0;

    while (std::getline(file, line) && row_idx < rows) {
        // '[' と ']' をスペースに置換する
        std::replace(line.begin(), line.end(), '[', ' ');
        std::replace(line.begin(), line.end(), ']', ' ');

        std::stringstream ss(line);
        std::string temp_val;
        int col_idx = 0;

        // 文字列として数値を取り出し、Scalar 型に変換
        while (ss >> temp_val && col_idx < cols) {
            std::stringstream val_ss(temp_val);
            Scalar val;
            val_ss >> val;
            
            B(row_idx, col_idx) = val;
            col_idx++;
        }
        row_idx++;
    }

    return B;
}