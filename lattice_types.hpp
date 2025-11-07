#pragma once

#include <Eigen/Dense>

// DeepLLLで使用する ScalarType を long double に統一
using ScalarType = long double; 
using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;