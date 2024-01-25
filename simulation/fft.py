from math import pi
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
isFit = False
isShow = False
FitByHand = False
fit_range_by_hand = [0.e12,.3e12]
fit_height = 30
fit_width = .5e11
fit_range_width = .8e11
plot_range = [0, 1*10**12]
plot_range_wide = [0., 2*10**12]

# data labels
PHASE_CONTRIB = "mix"
MAIN_SUB = "mix"
ANGLE_DELTA_DEG = "30"
TIME_MIN = "10"
ANGLE_FROM_PARALLEL_DEG = "5e-2"
LMD_USED_MIN = "7e-10"
LMD_USED_MAX = "10e-10"

output_extension="png"

# data files
FORMAT = PHASE_CONTRIB + "_" + MAIN_SUB + "_" + ANGLE_DELTA_DEG + \
    "deg_" + TIME_MIN + "min_ALPHA" + ANGLE_FROM_PARALLEL_DEG + \
    "_lmd" + LMD_USED_MIN + "to" + LMD_USED_MAX
inputFile = "./dat/theoretical/" + FORMAT + ".dat"
outputFile = "./oscil_graph/theoretical/fourier/" + FORMAT + "."+output_extension
outputFileWide = "./oscil_graph/theoretical/fourier_wide/" + FORMAT + "."+output_extension

# fitting range calculation
g = 9.8
h = 6.62607015e-34
m = 1.675e-27
gap = 189e-6
mirror_distance = 150e-3
theta = 1.05 * pi / 180
delta = float(ANGLE_DELTA_DEG) * pi / 180
peak_position = 2 * pi * g * pow(m / h, 2) * 2 * gap * \
    mirror_distance / np.tan(2 * theta) * np.sin(delta)
# (O-H)/(O+H)


def eq(x, y): return (x-y)/(x+y)

# fitting function


def gaussian(x, *params):
    return params[0] * np.exp(-(x - params[1]) ** 2 / (2 * params[2] ** 2)) + params[-1]


# Load data from file
inputList = np.loadtxt(inputFile, delimiter=' ')

wavelength = list(inputList[:, 0])
data = list(map(eq, inputList[:, 1], inputList[:, 2]))

# Fourier trans
amp = np.abs(np.fft.fft(data))
freq = 2*pi*np.abs(np.fft.fftfreq(len(wavelength),
                   np.mean(np.diff(wavelength))))

# plot
plt.figure(figsize=(8, 5))
plt.title('Fourier Transform ' + PHASE_CONTRIB + " " + MAIN_SUB + ' at delta = ' + ANGLE_DELTA_DEG + " deg ALPHA" + ANGLE_FROM_PARALLEL_DEG)
plt.xlabel('wave number [/m]')
plt.ylabel('Magnitude')
plt.axvline(x=peak_position, color="red")
plt.text(peak_position, 0, "expected peak = {:.4e}".format(peak_position), color="red")
plt.plot(freq, amp)

# Fitting
if isFit:
    if FitByHand:
        fit_range = fit_range_by_hand
        peak_position = np.mean(fit_range)
        fit_width = fit_range_by_hand[1] - fit_range_by_hand[0]
    else:
        fit_range = [peak_position - fit_range_width, peak_position + fit_range_width]
    print("peak position is {:e}".format(peak_position))
    freq_fit = freq[(freq > fit_range[0]) & (freq < fit_range[1])]
    amp_fit = amp[(freq > fit_range[0]) & (freq < fit_range[1])]
    popt, pocv = curve_fit(gaussian, freq_fit, amp_fit, p0=[
                           fit_height, peak_position, fit_width, 0])
    perr = np.sqrt(pocv)
    print('peak at {:e}'.format(popt[1]))
    fit = gaussian(freq_fit, *popt)
    plt.plot(freq_fit, fit)
    plt.text(popt[1]+fit_width, popt[0]*0.8,
             "highest freq = {:.4e} +- {:.4e}".format(popt[1], popt[2]))

plt.xlim(plot_range[0], plot_range[1])
# plt.ylim(-200, 2*10**3)
plt.grid(False)
plt.savefig(outputFile)
plt.close()

plt.title('Fourier Transform ' + PHASE_CONTRIB + " " + MAIN_SUB + ' at delta = ' + ANGLE_DELTA_DEG + " deg ALPHA" + ANGLE_FROM_PARALLEL_DEG)
plt.xlabel('wave number [/m]')
plt.ylabel('Magnitude')

plt.xlim(plot_range_wide[0], plot_range_wide[1])
plt.ylim(-1, 10)
plt.plot(freq, amp)
plt.savefig(outputFileWide)

if isFit:
    plt.plot(freq_fit, fit)
    plt.text(popt[1]+fit_width, popt[0]*0.8,
             "highest freq = {:e} +- {:e}".format(popt[1], popt[2]))

if isShow:
    plt.show()
# plt.close()
