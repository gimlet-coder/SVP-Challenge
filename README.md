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
    * DeepLLL / MLLL アルゴリズム
    * BKZ アルゴリズム
* **C++による数値計算・アルゴリズム実装能力の向上**
* **Git/GitHubを用いたバージョン管理の実践**

---

## 3. 使用技術

* **言語:** C++ (C++17)
* **コンパイラ:** g++ (MSYS2 UCRT64 環境推奨)
* **ライブラリ:** **Eigen** (行列演算), **Boost** (多倍長演算)

---

## 4. 実装済み機能

* `Gram_Schmidt.cpp`: グラム・シュミットの直交化
* `lagrange_basis_reduction.cpp`: ラグランジュの基底簡約アルゴリズム
* `Size-reduce.cpp`: サイズ基底簡約アルゴリズム
* `LLL.cpp`: LLL基底簡約アルゴリズム 
* `DeepLLL.cpp`: DeepLLL基底簡約アルゴリズム
* `MLLL.cpp`: MLLL基底簡約アルゴリズム　(β version 一次従属の行列に非対応)

**Upgrade**: 多倍長浮動小数対応
    → 高次元･高精度の格子基底簡約を正確に実行可能

---
## 5. 実行方法

本プロジェクトは MSYS2 UCRT64 環境でのビルドを前提としています。

1.  **依存ライブラリのインストール**
    
    格子計算の精度向上のため、Eigen に加えて Boost をインストールする必要があります。

    MSYS2の 「UCRT64」 ターミナルを開き、以下のコマンドでライブラリをインストールしてください。 
    (※従来のMinGW64版とはパッケージ名が異なりますのでご注意ください)
    
    ```bash
    # Eigen (行列演算), Boost (ユーティリティ) をインストール
    pacman -S mingw-w64-ucrt-x86_64-eigen3 mingw-w64-ucrt-x86_64-boost
    # GMP (任意精度整数) と MPFR (任意精度浮動小数点数) をインストール
    pacman -S mingw-w64-ucrt-x86_64-gmp mingw-w64-ucrt-x86_64-mpfr
    ```

2.  **ライブラリパスの設定(必要に応じて)**
    Makefile はデフォルトで C:/msys64/ucrt64/include/... を参照するように設定されています。 MSYS2を標準的な場所にインストールしている場合、編集は不要です。

    もし異なる場所にインストールしている場合は、Makefile 内の以下の行を環境に合わせて書き換えてください。

    ```Makefile
    EIGEN_PATH = C:/(各々の環境に応じて編集)/eigen3
    BOOST_PATH = C:/(各々の環境に応じて編集)/include
    ```

3.  **リポジトリのクローン**
    
    ```bash
    git clone [https://github.com/gimlet-coder/SVP-Challenge.git](https://github.com/gimlet-coder/SVP-Challenge.git)
    cd SVP-Challenge
    ```

4.  **コンパイル**
    ターミナル (PowerShell または MSYS2) で以下のコマンドを実行します。

    ```bash
    # コンパイル (Makefileが自動的に処理します)
    make

    # 'make' コマンドが見つからないと言われた場合
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

* [ ] MLLLアルゴリズムの完全な実装 (←Next)
* [ ] 格子上の最短ベクトルの数え上げ (ENUM) の実装
* [ ] BKZアルゴリズムの実装
* [ ] LLLアルゴリズムの改良
* [ ] ベクトル計算の効率化
* [ ] 100次元程度のSVPチャレンジに挑戦
* [ ] プログラムの並列化実装