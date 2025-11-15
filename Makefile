# コンパイラを指定
CXX = g++

# コンパイルフラグ: C++17, 警告をすべて表示, デバッグ情報付き
CXXFLAGS = -std=c++17 -Wall -g

# --- ★ここを自分の環境に合わせて編集★ ---
# Eigenライブラリのヘッダがある "親" フォルダへのパス
EIGEN_PATH = C:/msys64/mingw64/include/eigen3


# Boostライブラリのヘッダがある "親" フォルダへのパス
BOOST_PATH = C:/msys64/mingw64/include
# --- ★編集ここまで★ ---

# インクルードパスの指定
INCLUDES = -I$(EIGEN_PATH) -I$(BOOST_PATH)

# 実行可能ファイルの名前
TARGET = svp_solver

# --- ★ここを自分のファイル構成に合わせて編集★ ---
# コンパイル対象の .cpp ファイルをすべて列挙
# (例: main.cpp を追加し、Size-reduce_partial.cpp も追加)
SRCS =main.cpp Gram_Schmidt.cpp lagrange_basis_reduction.cpp Size-reduce.cpp Size-reduce_partial.cpp GSOUpdate-LLL_partial.cpp LLL.cpp GSOUpdate-DeepLLL_partial.cpp DeepLLL.cpp
# --- ★編集ここまで★ ---

# .cpp から .o を生成する
OBJS = $(SRCS:.cpp=.o)

# デフォルトのターゲット ( `make` とだけ打った時に実行される)
# ターゲットに .exe を明記する
all: $(TARGET).exe

# 実行可能ファイルを生成するルール (ターゲットに .exe を明記)
$(TARGET).exe: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET).exe $(OBJS)

# .cpp から .o (オブジェクトファイル) を生成するルール
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# `make clean` で実行されるルール (Windowsの 'del' コマンドに変更)
clean:
	del /F $(TARGET).exe $(OBJS)