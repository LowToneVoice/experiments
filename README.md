# 2023 P2 課題研究 COW実験

## 現状のおおまかな課題

詳細はissuesに記載

### 数値計算

- [x] ROOT Cern上でのFFT実装
  - [x] tutorialの確認
    - [x] tutorialのFFTがFFTされたように見えない
  - [x] 実装
- [x] JPARCのROOTファイルに沿った実装
  - [x] positionの実装
  - [x] lambda_roi_histogramの実装
  - [x] oscillation_lambdaの変更
  - [x] read_lambdaの変更
  - [x] sim_lambdaの変更
- [ ] 誤差の実装
  - [x] lambdaの誤差の組入
  - [ ] lambdaの誤差の確認
  - [ ] アラインメントの誤差の組入
  - [ ] アラインメントの誤差の確認
  - [x] 統計誤差の組入
  - [ ] 統計誤差の確認
- [x] 統計数の把握

### 実験装置

- [ ] 設計
  - [x] 要求の整理
  - [ ] 図面作成
- [ ] 必要機材の整理
  - [x] エタロン
  - [x] ステージ
  - [ ] 試料とステージのインターフェース
  - [ ] アラインメント配置の測定機器
    - [ ] 長さ計測
    - [ ] 角度計測
- [ ] 必要機材の調達確認
  - [ ] エタロン
  - [ ] ステージ
  - [ ] エタロンとステージのインターフェース
- [ ] 必要機材の稼働確認
  - [ ] ステージ操作方法
- [ ] 組み立て

## 実験概要

基本的にintroduction下「やりたい実験.pdf」に記載。
うまく実験が回ったらtexに。

## フォルダ構造・命名規則

``` folder tree
.
├── introduction
│   ├── やりたい実験.pdf: 実験ノートみたいな使い方してる
│   └── presentation.pptx: 実験概要プレゼン
│
├── papers: 参考文献 (gitには上げない)
│
├── phase-calc: 位相の概算に使用
│   ├── Data_all_mikata.md: Data_all.datのフォーマット及び入射角1.05ºでの位相 (誤差付き)
│   ├── phase.cpp: 位相とその誤差の概算
│   ├── prob_error.gpl: 波動関数から出る波長-確率分布の誤差を概算
│   ├── prob_error.pdf: 波動関数から出る波長-確率分布の誤差
│   └── ...
│
├── simulation
│   ├── beam_count: 波長-ビーム強度ヒストグラム
│   │   ├── montecarlo: モンテカルロシミュレーション
│   │   │   ├── normal: 全波長
│   │   │   └── zoom: 注目領域のみ
│   │   ├── position: TDC2次元プロット
│   │   └── theoretical: 波動関数から求めたO/H-beam観測確率
│   ├── BL05: ビームラインの波長-強度関係 (J-PARCより)
│   ├── chisq: グラフの振動フィッティングの精度計算 (不使用)
│   ├── dat: シミュレーションの計算結果生データ
│   │   ├── experiments: 実験データ
│   │   ├── montecarlo: モンテカルロシミュレーション
│   │   │   ├── dat: テキストファイル
│   │   │   └── root: rootファイル
│   │   └── theoretical: 波動関数から直接計算した確率分布
│   ├── oscil_graph: (I_H-I_O)/(I_H+I_O) シミュレーション結果のグラフ
│   │   ├── montecarlo: モンテカルロシミュレーション
│   │   │   ├── fourier: 振動のフーリエ変換
│   │   │   ├── fourier_wide: 振動のフーリエ変換 (波長範囲最大)
│   │   │   ├── normal: 振動グラフ
│   │   │   └── zoom: 振動グラフの拡大
│   │   └── theoretical: 波動関数から直接計算した確率分布
│   ├── reflection: ミラー反射率
│   │   ├── 2021P2: 2021年P2のコードと結果
│   │   ├── dat: シミュレーション結果テキストファイル
│   │   ├── graphs: シミュレーション結果
│   │   │
│   │   ├── reflection.cpp: 波長-反射率の1次元グラフ
│   │   └── reflection2d.cpp: 波長・入射角-反射率の2次元グラフ
│   ├── test: 挙動確認用。実験に直接は使わない
│   │
│   ├── fft.py: theoretical datファイルをFourier trans.してそれっぽいところをガウスフィッティング
│   ├── lambda_roi_histogram.cpp: 2次元位置を指定して波長-カウントヒストグラムを作成
│   ├── main.cpp: とりあえず動かすコード
│   ├── montecarlo.h: モンテカルロ計算に
│   ├── position.cpp: ビーム2次元位置のカウント数をカラープロット
│   ├── RESULT.md: 結果について
│   ├── simulate.h: シミュレーションに
│   └── ...
│
├── trash: ゴミ箱 (gitにはあげない)
│
├── .gitignore: git に追加しないファイルを指定
├── README.md: これ
└── ...
```

### 解析TTree

| Branch | 内容 |
| - | - |
| channel | H beam: 0, O beam: 1 |
| lambda | tofから換算した波長中心値 |

## 解析の流れ

### シミュレーション

1. `main.cpp` でシミュレーション
2. `fft.py` でフーリエ変換の理論予測

### 本実験

1. 実験
2. `position.cpp` でビームの2次元位置を確認
   1. 入出力ファイル変更 (l. 7-)
   2. 実行
3. `lambda_roi_histogram.cpp` で2次元位置から各ビームの波長-カウントヒストグラムを作成
   1. 入出力ファイル変更 (l. 12-)
   2. 切り取り位置の指定 (l. 23-)
   3. 実行
4. `` で振動を描画
5. `` で振動をFFT
   1. この中で重力加速度まで換算したい

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

- beam line: BL05 低発散ビームブランチ [J-PARC MLF](https://mlfinfo.jp/ja/bl05/)
- beam height 30e-3 (cf. 2024/01/19 thread)
- beam width 6e-3 (cf. 2024/01/19 thread)
- beam intensity 5.4e4 n/s/cm2
- ethalone area 12e-3 x 12e-3
- double slit width (sum of 2) 260e-6
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
- slit width 0.2e-3
- slit height 10e-3
- 縦方向発散 10 mrad (full open)

### チャンネル

- TDC H beam: 0
- TDC O beam: 1
- ADC H beam: 2
- ADC O beam: 3

## コードの動かし方

### Eigen/Denseを含むコードの動かし方

`g++ -g -o out test.cpp -std=c++14 -I /usr/local/include/eigen3`

ROOT も同時に動かす場合は
`g++ -o out main.cpp -std=c++17 -I /usr/local/include/eigen3 -I /usr/local/include/root -L /usr/local/lib/root -lCore -lImt -lRIO -lNet -lTree -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lThread -pthread`
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
