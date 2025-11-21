# コンパイラ (単に g++ とすることで、パスの通ったコンパイラを使用する)
CXX = g++

# コンパイルフラグ
CXXFLAGS = -std=c++17 -Wall -g

# --- ライブラリのパス設定 (デフォルト値) ---
# ?= を使うことで、ユーザーが make EIGEN_PATH=... と指定した場合に上書き可能になります
# MSYS2 UCRT64/MINGW64 の標準的なインストール先を指定

# Eigen3: 通常は include/eigen3 
EIGEN_PATH ?= C:/msys64/ucrt64/include/eigen3

# Boost: 通常は include 
BOOST_PATH ?= C:/msys64/ucrt64/include

# Linux 環境でも問題なく動作するように対応予定


# インクルードパス

INCLUDES = -I$(EIGEN_PATH) -I$(BOOST_PATH)

# 実行可能ファイル名
TARGET = svp_solver

# ソースファイル
SRCS = main.cpp Gram_Schmidt.cpp lagrange_basis_reduction.cpp Size-reduce.cpp Size-reduce_partial.cpp GSOUpdate-LLL_partial.cpp LLL.cpp GSOUpdate-DeepLLL_partial.cpp DeepLLL.cpp MLLL.cpp ENUM.cpp

# オブジェクトファイル
OBJS = $(SRCS:.cpp=.o)

# ターゲット
all: $(TARGET).exe

$(TARGET).exe: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET).exe $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	del /F $(TARGET).exe $(OBJS)