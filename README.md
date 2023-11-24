# 2023 P2 課題研究 COW実験

## 実験概要

### 背景

物質波の時間発展は Schrödinger eq. に従うが、非相対論的にはハミルトニアンに重力ポテンシャル $mgx$ が加わる。
重力ポテンシャルが異なる2つの経路を通ると位相差 $-2\pi\lambda m^2gA/h^2$ が生じる。
この効果は Overhauser と Colella により提唱され、Colella, Overhauser, Werner により実験的に検証された。

1988 年に Werner, Kaiser, Arif, Clothier により追試、まとめられたところによると、


## ファイル構造

## 定数表

### 数学・物理定数

原則SI単位系でコーディングする

- pi 3.14159265359
- Planck const. 6.62607015e-34
- neutron mass 1.6749e-27
- elementary charge 1.60218e-19
- grav. acc. 9.8
- Optical pot. Ni 224.e-9 eV
- Optical pot. Ni -40.e-9 eV

### 装置設計の不可変定数

- DAQ freq. 62.5e+6
- gap thickness 189e-6
- wavelength min. 2e-10
- wavelength max. 9e-10
- layer thickness of Ni 13.35e-9
- layer thickness of Ti 9.83e-9
- bilayer count 8

### 装置設計の可変変数

- angle 1.05º
- beam time 1 h
- DAQ downsizing 16

## ROOT Cernの動かし方

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
- [Werner, Kaiser, Arif, Clothier, Physica B+C, Volume 151, Issues 1–2, 1988.](https://doi.org/10.1016/0378-4363(88)90141-6) (COWの新しい結果)
- [Y. Seki, 2011](http://hdl.handle.net/2433/142371) (中性子干渉実験周辺、特に屈折率計算)

### 使わないやつ

- [Bonse and Wroblewski, Phys. Rev. Lett. 51, 1401 – Published 17 October 1983](https://doi.org/10.1103/PhysRevLett.51.1401) (非慣性系における中性子干渉測定)
- [Bonse and Wroblewski, Phys. Rev. D 30, 1214 – Published 15 September 1984](https://doi.org/10.1103/PhysRevD.30.1214) (同)
- [arXiv:1701.00259](https://arxiv.org/abs/1701.00259) (Galiautdinov and Ryder, 2017. 相対論的COW)
- [S. A. Werner, J. -L. Staudenmann, and R. Colella Phys. Rev. Lett. 42, 1103 – Published 23 April 1979](https://doi.org/10.1103/PhysRevLett.42.1103) (Sagnac effect)
