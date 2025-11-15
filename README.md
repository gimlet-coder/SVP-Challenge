# C++によるSVPチャレンジ実装

[![Textbook](https://img.shields.io/badge/Reference-格子暗号解読-blue.svg)](https://amzn.asia/d/40cFynK)

C++を用いて、教科書『[格子暗号解読 ―SVPチャレンジ・LWEチャレンジ ―](https://amzn.asia/d/40cFynK)』を参考にしながら最短ベクトル問題(SVP)チャレンジに取り組むリポジトリです。

---

## 1. 概要

格子暗号解読の基礎となるSVP（Shortest Vector Problem）の求解アルゴリズムをC++で実装することを目標としています。



---

## 2. プロジェクトの目的

* **SVPおよび関連アルゴリズムの理解:**
    * グラム・シュミットの直交化
    * LLL (Lenstra–Lenstra–Lovász) アルゴリズム
* **C++による数値計算・アルゴリズム実装能力の向上**
* **Git/GitHubを用いたバージョン管理の実践**

---

## 3. 使用技術

* **言語:** C++ (C++17)
* **コンパイラ:** g++
* **ライブラリ:** **Eigen** (行列演算), **Boost** (ユーティリティ), **GMP/MPFR** (多倍長浮動小数点数演算)

---

## 4. 実装済み機能

* `Gram_Schmidt.cpp`: グラム・シュミットの直交化
* `lagrange_basis_reduction.cpp`: ラグランジュの基底簡約アルゴリズム
* `Size-reduce.cpp`: サイズ基底簡約アルゴリズム
* `LLL.cpp`: LLL基底簡約アルゴリズム 
* `DeepLLL.cpp`: DeepLLL基底簡約アルゴリズム

**Upgrade** 多倍長浮動小数対応
    → 高次元･高精度の格子基底簡約を正確に実行可能

---
## 5. 実行方法

本プロジェクトのコンパイルには以下ライブラリが必要です。

1.  **依存ライブラリのインストール**
    
    格子計算の精度向上のため、Eigen、Boost、GMP および MPFR が必須です。

    WindowsのMinGW環境でビルドする場合、MSYS2 MinGW x64 シェルで以下のコマンドを実行してインストールしてください。
    
    ```bash
    # Eigen (行列演算), Boost (ユーティリティ) をインストール
    pacman -S mingw-w64-x86_64-eigen3 mingw-w64-x86_64-boost
    # GMP (任意精度整数) と MPFR (任意精度浮動小数点数) をインストール
    pacman -S mingw-w64-x86_64-gmp mingw-w64-x86_64-mpfr
    ```

2.  **ライブラリパスの設定(必須)**
    プロジェクトルート（SVP-Challenge フォルダ直下）にある Makefile をテキストエディタで開いてください。

    EIGEN_PATH = および BOOST_PATH = で始まる行を見つけ、= の右側を、あなたのPC環境に合わせて必ず書き換えてください。

    (設定例: EIGEN_PATH = C:/msys64/mingw64/include/eigen3, BOOST_PATH = C:/msys64/mingw64/include)

3.  **リポジトリのクローン**
    
    ```bash
    git clone [https://github.com/gimlet-coder/SVP-Challenge.git](https://github.com/gimlet-coder/SVP-Challenge.git)
    cd SVP-Challenge
    ```

4.  **コンパイル**
    ターミナルで以下のコマンドを実行します。Makefile がビルドプロセスを自動化します。

    (WindowsのMinGW環境では mingw32-make が正式なコマンド名の場合があります)
    ```bash
    # (PowerShellで 'make' のエイリアス設定が完了している場合)
    make

    # 'make' が見つからないと言われた場合
    mingw32-make
    ```

5.  **実行**
    コンパイルが成功すると、フォルダ内に svp_solver.exe が生成されます。 以下のコマンドで実行してください。
    ```bash
    ./svp_solver.exe
    ```
---

6. **クリーンアップ(任意)**
    生成された実行ファイル (.exe) やオブジェクトファイル (.o) を削除する場合は、以下のコマンドを実行します。
    ```bash
    make clean
    # (または mingw32-make clean)
    ```

## 7. 更新予定

* [ ] MLLL基底簡約アルゴリズムの実装 (←Next)
* [ ] LLLアルゴリズムの改良
* [ ] ベクトル計算の効率化
* [ ] 100次元程度のSVPチャレンジに挑戦
* [ ] プログラムの並列化実装