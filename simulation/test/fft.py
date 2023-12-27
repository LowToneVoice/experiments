import numpy as np  # numpyをインポート
import matplotlib.pyplot as plt  # matplotlibをインポート

### 量子化(離散化) ###
f_s = 44100  # サンプリングレート f_s[Hz] (任意)
t_fin = 1  # 収録終了時刻 [s] (任意)
dt = 1/f_s  # サンプリング周期 dt[s]
N = int(f_s * t_fin)  # サンプル数 [個]

### 入力信号 y(t) ###
f0 = 440  # 周波数 f0[Hz]
t = np.arange(0, t_fin, dt)  # 離散時間 t[s]
y = np.sin(2*np.pi*f0*t)  # y(t) を生成

### y(t) のグラフ ###
plt.plot(t, y)  # y-t グラフのプロット
plt.xlim(0, 10*10**-3)  # 横軸に関する描画範囲指定
plt.show()  # グラフの表示

### FFT: tの関数をfの関数にする ###
y_fft = np.fft.fft(y)  # 離散フーリエ変換
freq = np.fft.fftfreq(N, d=dt)  # 周波数を割り当てる（※後述）
Amp = abs(y_fft/(N/2))  # 音の大きさ（振幅の大きさ）

### 音波のスペクトル ###
plt.plot(freq[1:int(N/2)], Amp[1:int(N/2)])  # A-f グラフのプロット
plt.xscale("log")  # 横軸を対数軸にセット
plt.show()
