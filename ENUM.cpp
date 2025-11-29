#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <limits>


#include "lattice_types.hpp"

bool ENUM(const Matrix &U, const Vector &B_norm, Scalar &R_square, Vector &v_out, const int k_begin, const int k_end, long long &node_count){
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

    while(true){ // step.9
        node_count++;
        const int global_k = k + k_begin - 1;
        Scalar diff = static_cast<Scalar>(v_temp(k - 1)) - c(k - 1);

        rho(k) = rho(k + 1) + diff * diff * B_norm(global_k); // step.10


        if(rho(k) <= R_square){ // step.11
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
                    R_square = rho(k);
                    FLAG_found = true;
std::cerr << "見つかったノルム: " << R_square << " (Nodes: " << node_count << ")" << std::endl;
                    if (v_temp(k - 1) >= c(k - 1)){
                        v_temp(k - 1) += w_int(k - 1);
                    }else{
                        v_temp(k - 1) -= w_int(k - 1);
                    }
                    w_int(k - 1)++;
                    continue;
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
