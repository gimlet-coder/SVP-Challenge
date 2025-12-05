# CMakeの最低限必要なバージョンを指定
cmake_minimum_required(VERSION 3.15)

# プロジェクト名を定義 
project(SVP_Challenge_DeepLLL)

# C++17 標準を使用するように設定
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------------------------------------------
# --- 実行ファイルの作成
# -------------------------------------------------------------
set(MAIN_SOURCES
    main.cpp
    LLL.cpp
    DeepLLL.cpp
    Gram_Schmidt.cpp
    GSOUpdate-DeepLLL_partial.cpp
    GSOUpdate-LLL_partial.cpp
    Size-reduce_partial.cpp
    MLLL.cpp
    ENUM.cpp
    BKZ.cpp
)
add_executable(main ${MAIN_SOURCES})
# -------------------------------------------------------------


# -------------------------------------------------------------
# --- 必要なライブラリを探す
# -------------------------------------------------------------
# Eigenの検索
find_package(Eigen3 CONFIG REQUIRED)

# Boostの検索 (multiprecisionモジュールのみに戻す)
find_package(Boost CONFIG REQUIRED COMPONENTS multiprecision)

# -------------------------------------------------------------


# -------------------------------------------------------------
# --- ライブラリのリンクとインクルード設定
# -------------------------------------------------------------
target_link_libraries(main PRIVATE 
    Eigen3::Eigen 
    Boost::multiprecision
)

# 【重要】GMP/MPFR/MPC のヘッダーパスの強制追加 (vcpkgのバグ回避策)
# MPFRのヘッダーエラーを回避するため、この行は残します
target_include_directories(main PRIVATE 
    ${Boost_INCLUDE_DIRS} #Boost のインクルードディレクトリを追加
    C:/vcpkg/installed/x64-windows/include
)

# (おまけ) Eigen が最適なCPU命令を使えるようにする
target_compile_features(main PRIVATE cxx_std_17)
target_compile_options(main PRIVATE /arch:AVX2) 
# -------------------------------------------------------------