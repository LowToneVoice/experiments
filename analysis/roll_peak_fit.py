import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ガウス関数の定義


def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2)


# データの手動設定
x_data = np.array([19,19.5,20,20.5])  # 手動でx座標の配列を調整
y_data = np.array([1,1,1,1])  # 手動でy座標の配列を調整

# ガウスフィッティングの実行
params, covariance = curve_fit(gaussian, x_data, y_data, p0=[1, 5, 1])
amplitude, mean, stddev = params

# プロット
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_data, gaussian(x_data, amplitude,
         mean, stddev), color='red', label='Fit')
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
