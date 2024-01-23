# ビームカウント見積

## 各種設定数値

### ビームライン情報

BL05 低発散ビームブランチを使用
[J-PARC MLF](https://mlfinfo.jp/ja/bl05/)

- beam height 30 mm
- beam width 6 mm (after collimation), 8 mm (before collimation)
- beam intensity 5.4e4 n/s/cm2
- beam発散角 5.4e-2 µstr
- 中性子ビーム成形B4Cスリット(7.5 m位置と12 m 位置)
- mix of wavelength $0.01\times(\lambda_\mathrm{max}-\lambda_\mathrm{min})$

### ビームカウントの波長依存性

bl05_namadata.ccでは

- mass = 1.675e-27
- hbar = 1.055e-34
- L_tof = 18022e-03

によって $\lambda=2.1959\times10^{-8}\tau$ となっている。
tofは 0 ~ 45000 µs を 1000 分割しているので、bin幅は
45000e-6 / 1000 * 2.1959e-8 = 9.88e-13 m.

| lambda (m) | count per second |
|-|-|
| 0.8e-9 | 0.22 |
| 0.6e-9 | 0.63 |

### エタロン情報

- ethalone height 12 mm
- ethalone width 12 mm
- double slit width (sum of 2) 260 µm
- layer thickness of Ni 13.35 nm (cf. Seki 2011, Tbl. 4.3)
- layer thickness of Ti 9.83 nm (ibid.)
- bilayer count 8

beam collimation外さないと30ºでエタロンがビームからはみ出る。
外せば収まる。

### 装置設計変数

- angle theta 1.05º
- angle delta 30º

## 計算

### beam area

(slit width) x (ethalone height) = 3.12e-2 cm2

## beam count at lambda region of interest

`d_lambda = 9.88e-13`,
then line in lambda-count histogram is `count_per_second = -2.05e9 * lambda + 1.86`.
Therefore, when using d_lambda,
`count_per_second = -2.05e9 * (d_lambda / 9.88e-13) * lambda + 1.86`.
