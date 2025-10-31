# C++によるSVPチャレンジ実装

[![Textbook](https://img.shields.io/badge/Reference-格子暗号解読-blue.svg)](https://amzn.asia/d/40cFynK)

C++を用いて、教科書『[格子暗号解読 ―SVPチャレンジ・LWEチャレンジ ―](https://amzn.asia/d/40cFynK)』を参考にしながら最短ベクトル問題(SVP)チャレンジに取り組むリポジトリです。

---

## 1. 概要

このプロジェクトは、格子暗号解読の基礎となるSVP（Shortest Vector Problem）の求解アルゴリズムをC++で実装し、その理解を深めることを目的としています。

学習と実装の過程を記録することで、就職活動用の技術ポートフォリオとして活用することも目指しています。

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
* **ライブラリ:** (Eigen/Dense)

---

## 4. 実装済み機能

* `Gram_Schmidt.cpp`: グラム・シュミットの直交化
* `lagrange_basis_reduction.cpp`: ラグランジュの基底簡約アルゴリズム
* `Size-reduce.cpp`: サイズ基底簡約アルゴリズム

* (今後実装する機能...)

---

## 5. 実行方法

1.  リポジトリをクローンします。
    ```bash
    git clone [https://github.com/gimlet-coder/SVP-Challenge.git](https://github.com/gimlet-coder/SVP-Challenge.git)
    cd SVP-Challenge
    ```

2.  (ここにコンパイル方法を具体的に書いてください)
    
    **コンパイル例:**
    ```bash
    # (例) g++ を使って全ての .cpp ファイルをコンパイルする場合
    g++ -std=c++17 -o svp_solver *.cpp
    ```

3.  実行します。
    ```bash
    ./svp_solver
    ```

---

## 6. 今後の課題・TODO

* [ ] LLLアルゴリズムの実装
* [ ] ベクトル計算の効率化
* [ ] 実際にSVPチャレンジに取り組む