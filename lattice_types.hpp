#pragma once

#include <boost/config.hpp> 
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/detail/default_ops.hpp>


#include <Eigen/Dense>




// --- 1. 多倍長精度 (Scalar) の型を定義 ---
using Scalar = boost::multiprecision::cpp_dec_float_100; // 100桁精度の多倍長浮動小数点数


const Scalar MULTI_PRECISION_EPSILON = Scalar("1e-90"); // 90桁までを信頼する

// --- 2. Eigen が Scalar 型を認識できるようにする「おまじない」 ---
namespace Eigen {
    template<> struct NumTraits<Scalar> : NumTraits<long double> {
        typedef Scalar Real;
        typedef Scalar NonInteger;
        typedef Scalar Nested;
        typedef Scalar Literal;
        typedef Scalar FloatingPoint;

        enum {
            IsComplex = 0, IsInteger = 0, IsSigned = 1, RequireInitialization = 1,
            ReadCost = 1, AddCost = 4, MulCost = 8
        };
    };
}

// --- 3. Matrix と Vector の型定義 ---
// 以前の long double を Scalar に置き換え
using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>; // 列ベクトル
using Matrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;


// --- 4. Scalar 型における演算の下準備 ---

using boost::multiprecision::abs;
using boost::multiprecision::round;
using boost::multiprecision::floor;
using boost::multiprecision::ceil;
using boost::multiprecision::sqrt;
using boost::multiprecision::pow;
