# 2023 P2 課題研究 COW実験

## 現状の課題

### 数値計算

[ ] ROOT Cern上でのFFT実装
[ ] JPARCのROOTファイルに沿った実装

### 実験装置

## 実験概要

基本的にintroduction下「やりたい実験.pdf」に記載。
うまく実験が回ったらtexに。

## フォルダ構造

``` folder tree
.
├── introduction
│   ├── やりたい実験.pdf：実験ノートみたいな使い方してる
│   └── presentation.pptx：実験概要プレゼン
│
├── papers：参考文献 (gitには上げない)
│
├── phase-calc：位相の概算に使用
│   ├── Data_all_mikata.md：Data_all.datのフォーマット及び入射角1.05ºでの位相 (誤差付き)
│   ├── phase.cpp：位相とその誤差の概算
│   ├── prob_error.gpl：波動関数から出る波長-確率分布の誤差を概算
│   ├── prob_error.pdf：波動関数から出る波長-確率分布の誤差
│   └── ...
│
├── simulation
│   ├── beam_count：波長-ビーム強度ヒストグラム
│   │   ├── montecarlo：モンテカルロシミュレーション
│   │   └── theoretical：波動関数から求めたO/H-beam観測確率
│   │       ├── noref：反射・透過係数を無視
│   │       └── ref：反射・透過係数を入れてる
│   ├── BL05：ビームラインの波長-強度関係 (KEKより)
│   ├── chisq：グラフの振動フィッティングの精度計算
│   │   ├── chi2：振動フィッティングの初期値-chi2関係
│   │   ├── k-value：振動フィッティングの初期値-フィッティング結果
│   │   └── oscillation：振動フィッティング
│   ├── dat：シミュレーションの計算結果生データ
│   │   ├── **montecarlo：モンテカルロシミュレーション**
│   │   │   ├── **noref：反射透過の位相を考慮していない**
│   │   │   │   ├── **dim2：入射角依存性を考慮に入れる**
│   │   │   │   └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   │   │   └── **ref：反射透過の位相を考慮**
│   │   │       ├── **dim2：入射角依存性を考慮に入れる**
│   │   │       └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   │   └── **theoretical：波動関数から直接計算した確率分布**
│   │       ├── **noref：反射透過の位相を考慮していない**
│   │       │   ├── **dim2：入射角依存性を考慮に入れる**
│   │       │   └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   │       └── **ref：反射透過の位相を考慮**
│   │           ├── **dim2：入射角依存性を考慮に入れる**
│   │           └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   ├── oscil_graph：(I_H-I_O)/(I_H+I_O) シミュレーション結果のグラフ
│   │   ├── **montecarlo：モンテカルロシミュレーション**
│   │   │   ├── **noref：反射透過の位相を考慮していない**
│   │   │   │   ├── **dim2：入射角依存性を考慮に入れる**
│   │   │   │   └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   │   │   └── **ref：反射透過の位相を考慮**
│   │   │       ├── **dim2：入射角依存性を考慮に入れる**
│   │   │       └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   │   └── **theoretical：波動関数から直接計算した確率分布**
│   │       ├── **noref：反射透過の位相を考慮していない**
│   │       │   ├── **dim2：入射角依存性を考慮に入れる**
│   │       │   └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   │       └── **ref：反射透過の位相を考慮**
│   │           ├── **dim2：入射角依存性を考慮に入れる**
│   │           └── **lambda：波長特性のみを考え入射角は1.05ºに固定**
│   ├── reflection：ミラー反射率
│   │   ├── 2021P2：2021年P2のコードと結果
│   │   │
│   │   ├── reflNiTi-sim2.pdf：1枚 (bilayer×8) の波長-反射率
│   │   └── reflNiTi2d.pdf：1枚の波長・角度-反射率
│   ├── test：挙動確認用。実験に直接は使わない
│   │
│   ├── fft.py：theoretical datファイルをFourier trans.してそれっぽいところをガウスフィッティング
│   ├── oscillation_lambda.cpp：ファイルから $(I_H-I_O)/(I_H+I_O)$ を計算
│   ├── read_2d.cpp：シミュレーション dat ファイルから波長・角度-強度関係をヒストグラムにする
│   ├── read_lambda.cpp：シミュレーション dat ファイルから波長-強度関係をヒストグラムにする
│   ├── sim_2d.cpp：波長・角度-強度関係をモンテカルロシミュレーション＆波動関数計算
│   ├── sim_lambda.cpp：波長-強度関係をモンテカルロシミュレーション＆波動関数計算
│   ├── theoretical_lambda.gpl : 波動関数から求めた確率分布と振動の描画
│   ├── theoretical_2d.gpl : 波動関数から求めた角度依存確率分布と振動の描画
│   └── ...
│
├── trash：ゴミ箱 (gitにはあげない)
│
├── .gitignore：git に追加しないファイルを指定
├── README.md：これ
└── ...
```

## 解析の流れ

### シミュレーション

1. `simulation/sim_lambda.cpp` でデータを生成
   1. 載せたい位相を選択 (l. 200あたり)
   2. 出力ファイルを変更 (l. 10-)
   3. 適宜数値を変更 (関数宣言の前にまとめてる)
   4. 実行
2. `simulation/read_lambda.cpp` でビームのヒストグラムとその理論分布を作成
   1. 入出力ファイルを変更 (l. 10-)
   2. 実行
3. `theoretical_lambda.gpl` で振動の理論分布を描画
   1. 入出力ファイルを変更 (l. 2-3)
   2. 適宜範囲を変更 (l. 4-7)
   3. 実行
4. `oscillation_lambda.cpp` で振動のモンテカルロシミュレーションを実行
   1. 入出力ファイルを変更 (l. 10-)
   2. 実行

### 本実験

## 設定表

### 数学・物理定数

原則SI単位系でコーディングする

- pi 3.14159265359
- Planck const. 6.62607015e-34
- Avogadro const. 6.02214e23
- neutron mass 1.6749e-27
- elementary charge 1.60218e-19
- grav. acc. 9.8
- Optical pot. Ni 224.e-9 eV (Seki 2011, Tbl. 4.3)
- Optical pot. Ti -40.e-9 eV (ibid.)
- Optical pot. SiO2 90.5e-9 eV (ibid.)
- atomic mass of Ni 58.6934 g/mol
- atomic mass of Ti 47.867 g/mol
- mass density of Ni 8.908 g/cm3
- mass density of Ti 4.506 g/cm3
- neutron-nuclear scattering length of Ni 10.3 fm
- neutron-nuclear scattering length of Ti -3.36 fm
- atomic density of Si 5.00e22 (op. cit. p. 82)

### 装置設計の不可変定数

- DAQ freq. 62.5e+6
- gap thickness $(189\pm0.1)\times10^{-6}$
- wavelength min. 2e-10
- wavelength max. 9e-10
- layer thickness of Ni 13.35e-9 (op. cit. Tbl. 4.2)
- layer thickness of Ti 9.83e-9 (ibid.)
- bilayer count 8
- mix of wavelength $0.01\times(\lambda_\mathrm{max}-\lambda_\mathrm{min})$

### 装置設計の可変変数

- mirror distance 150e-3
- angle $(1.05\pm0.1)\degree$ (あまり変えたくない)
- beam time 1 h
- DAQ downsizing 16
- used wavelength min. 6.9e-10
- used wavelength min. 8.4e-10
- theta min. 0.2º
- theta max. 1.5º
- theta interval 0.01º
- total length $1.000\pm.005$
- delta 30º

### チャンネル

- TDC H beam: 0
- TDC O beam: 1
- ADC H beam: 2
- ADC O beam: 3

## コードの動かし方

### Eigen/Denseを含むコードの動かし方

`g++ -g -o out test.cpp -std=c++14 -I /usr/local/include/eigen3`

ROOT も同時に動かす場合は
`g++ -g -o simulation_code test.cpp -std=c++14 -I /usr/local/include/eigen3 -I /usr/local/include/root -L /usr/local/lib/root -lCore -lImt -lRIO -lNet -lTree -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lThread -pthread`
のようにする。

### TTree の対話型

``` ROOT Cern
$ root tree.root
...
root [1] _file0->ls()   // rootファイルの中身を一覧
TFile**         tree.root
 TFile*         tree.root
   KEY: TTree   tree;234        tree [current cycle]
   KEY: TTree   tree;233        tree [backup cycle]
root [2] tree->Print()  // 情報を表示
******************************************************************************
*Tree    :tree      : tree                                                   *
*Entries : 1105890000 : Total =     13271238886 bytes  File  Size = 6995942660 *
*        :          : Tree compression factor =   1.90                       *
******************************************************************************
*Br    0 :lambda    : lambda/D                                               *
*Entries :1105890000 : Total  Size= 8847425450 bytes  File Size  = 6975232995 *
*Baskets :     3281 : Basket Size=    4267008 bytes  Compression=   1.27     *
*............................................................................*
*Br    1 :channel   : channel/I                                              *
*Entries :1105890000 : Total  Size= 4423813079 bytes  File Size  = 20649327 *
*Baskets :     2689 : Basket Size=    2133504 bytes  Compression= 214.23     *
*............................................................................*
root [3] tree->Draw("lambda>>h1 (10000,2e-10,9e-10)","channel==0")  // ヒストグラムにしたいブランチを選択
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
```

## 参考文献

### 使うやつ

- [R. Colella, A. W. Overhauser, and S. A. Werner Phys. Rev. Lett. 34, 1472 – Published 9 June 1975](https://doi.org/10.1103/PhysRevLett.34.1472) (COW実験)
- [Fujie et. al. 2023. arXiv:2308.01922](https://doi.org/10.48550/arXiv.2308.01922) (多層膜ミラー)
- [A. W. Overhauser and R. Colella, Phys. Rev. Lett. 33, 1237 – Published 11 November 1974](https://doi.org/10.1103/PhysRevLett.33.1237) (COW theory)
- [多層膜ミラーを用いた中性子干渉計の作成と重力加速度の測定](https://www-he.scphys.kyoto-u.ac.jp/gakubu/P2/P2-21/P2_2021_report_neutron.pdf) (2021 P2 report)
- [多層膜ミラーを用いた中性子干渉計の作成と重力加速度の測定](https://www-he.scphys.kyoto-u.ac.jp/gakubu/P2/P2-21/P2_2021_slide_neutron.pdf) (2021 P2 slides)
  - cf. [反射率計算コード](https://drive.google.com/drive/folders/1OXl9TjSrukBzPKXic_n39EXW2YIPZwYu?usp=drive_link)
- [Werner, Kaiser, Arif, Clothier, Physica B+C, Volume 151, Issues 1–2, 1988.](https://doi.org/10.1016/0378-4363(88)90141-6) (COWの新しい結果)
- [Y. Seki, 2011](http://hdl.handle.net/2433/142371) (中性子干渉実験周辺、特に屈折率計算)

### 使わないやつ

- [Bonse and Wroblewski, Phys. Rev. Lett. 51, 1401 – Published 17 October 1983](https://doi.org/10.1103/PhysRevLett.51.1401) (非慣性系における中性子干渉測定)
- [Bonse and Wroblewski, Phys. Rev. D 30, 1214 – Published 15 September 1984](https://doi.org/10.1103/PhysRevD.30.1214) (同)
- [arXiv:1701.00259](https://arxiv.org/abs/1701.00259) (Galiautdinov and Ryder, 2017. 相対論的COW)
- [S. A. Werner, J. -L. Staudenmann, and R. Colella Phys. Rev. Lett. 42, 1103 – Published 23 April 1979](https://doi.org/10.1103/PhysRevLett.42.1103) (Sagnac effect)
