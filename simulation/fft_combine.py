from math import pi
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

plot_range = [0, 1*10**12]
plot_range_wide = [0., 2*10**12]

# data labels
ANGLE_DELTA_DEG = "0"
TIME_MIN = "10"
LMD_USED_MIN = "7e-10"
LMD_USED_MAX = "8e-10"
ALPHAS = ["0e-3", "-1e-3", "1e-3", "3e-3", "-3e-3"]
styles = ["-", ":", ":", ":", ":"]
grayCol = 0.

# fitting range calculation
g = 9.8
h = 6.62607015e-34
m = 1.675e-27
gap = 189e-6
mirror_distance = 150e-3
theta = 1.05 * pi / 180
delta = float(ANGLE_DELTA_DEG) * pi / 180
peak_position = 2*pi*g*pow(m/h, 2)*2*gap * \
    mirror_distance/np.tan(2*theta)*np.sin(delta)


# (H-O)/(H+O)
def eq(x, y): return (x-y)/(x+y)

# fitting fnc


def gaussian(x, *params):
    return params[0] * np.exp(-(x - params[1]) ** 2 / (2 * params[2] ** 2)) + params[-1]


# plot
plt.figure(figsize=(8, 5))
plt.title('Fourier Transform at delta = '+ANGLE_DELTA_DEG+' deg by alpha')
plt.xlabel('wave number [/m]')
plt.ylabel('Magnitude')
plt.axvline(x=peak_position, color="red")
plt.text(peak_position, 0, "expected peak = {:.4e}".format(
    peak_position), color="red")
plt.ylim(-1, 100)
plt.xlim(plot_range[0], plot_range[1])
plt.grid(False)

for (ANGLE_FROM_PARALLEL_DEG, sty) in zip(ALPHAS, styles):
    # data files
    FORMAT = "mix_mix_"+ANGLE_DELTA_DEG+"deg_"+TIME_MIN+"min_ALPHA" + \
        ANGLE_FROM_PARALLEL_DEG+"_lmd"+LMD_USED_MIN+"to"+LMD_USED_MAX
    inputFile = "./dat/theoretical/"+FORMAT+".dat"

    # load data from file
    inputList = np.loadtxt(inputFile, delimiter=' ')
    wavelength = list(inputList[:, 0])
    data = list(map(eq, inputList[:, 1], inputList[:, 2]))

    # Fourier trans
    amp = np.abs(np.fft.fft(data))
    freq = 2*pi*np.abs(np.fft.fftfreq(len(wavelength),
                       np.mean(np.diff(wavelength))))

    # plot
    plt.plot(freq, amp, sty, label=ANGLE_FROM_PARALLEL_DEG)
    # grayCol += .3

plt.legend(loc=1)
outputFile = "./oscil_graph/theoretical/fourier/delta"+ANGLE_DELTA_DEG+".pdf"
plt.savefig(outputFile)
plt.close()
