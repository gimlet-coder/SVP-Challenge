#pragma once

#include <boost/config.hpp> 
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp> //追加: 多倍長整数
#include <boost/multiprecision/detail/default_ops.hpp>


#include <Eigen/Dense>
#include <string>





// --- 1. 多倍長精度 (Scalar) の型を定義 ---
// 基底用 多倍長整数 (Arbitrary Precision Integer)
using Integer = boost::multiprecision::cpp_int;

// GSO/計算用 多倍長浮動小数点数 (Arbitrary Precision Float)
using FloatBackend = boost::multiprecision::cpp_dec_float<200>; // ここで有効数字を決める
using Real = boost::multiprecision::number<FloatBackend>; 

// 信頼区間 (Float用)
const Real MULTI_PRECISION_EPSILON = Real("1e-190");

// --- 2. Eigen が Scalar 型を認識できるようにする「おまじない」 ---
namespace Eigen {
    // Integer (cpp_int) 用の設定
    template<> struct NumTraits<Integer> : GenericNumTraits<Integer> {
        typedef Integer Real;// 整数の「実数部」も整数として扱う
        typedef Integer NonInteger;
        typedef Integer Nested;
        typedef Integer Literal;
        
        enum {
            IsComplex = 0,
            IsInteger = 1,// 整数であることを明示
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1, AddCost = 4, MulCost = 8
        };
    };

    // Real (cpp_dec_float) 用の設定
    template<> struct NumTraits<Real> : GenericNumTraits<Real> {
        typedef Real Real;
        typedef Real NonInteger;
        typedef Real Nested;
        typedef Real Literal;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1, AddCost = 4, MulCost = 8
        };
    };
}
// --- 3. Matrix と Vector の型定義 ---
// 以前の long double を Scalar に置き換え

using IntVector = Eigen::Matrix<Integer, Eigen::Dynamic, 1>;
using IntMatrix = Eigen::Matrix<Integer, Eigen::Dynamic, Eigen::Dynamic>; // 基底行列には IntMatrix を使う

using RealVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using RealMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>; // GSOとかの計算は RealMatrix を使う


// --- 4. Scalar 型における演算の下準備 ---

using boost::multiprecision::abs;
using boost::multiprecision::round;
using boost::multiprecision::floor;
using boost::multiprecision::ceil;
using boost::multiprecision::sqrt;
using boost::multiprecision::pow;


// --- 5. lattice file を読み取る関数 ---
IntMatrix load_challenge_matrix(const std::string &filename, int rows, int cols);
/**
 * @brief filename の lattice を読み込んで B として返す関数
 * @param filename [in] lattice の保存しているファイル名を入れる
 * @param rows, cols [in] それぞれ行数と列数を表す
 */