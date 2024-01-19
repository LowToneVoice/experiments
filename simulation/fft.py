from math import pi
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
isFit = False
plot_range = [0.1*10**11, .7*10**12]
fit_range = [0.2*10**12, .7*10**12]

inputFile = "./dat/theoretical/ref/lambda/mix_mix_90deg_N8e6_ALPHA1e-3.dat"
outputFile = "./oscil_graph/theoretical/ref/lambda/mix_mix_90deg_N8e6_ALPHA1e-3_fourier.pdf"

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
# print(len(freq))
# print(freq[50:100])
# print(amp)

plt.figure(figsize=(8, 5))
plt.title('Fourier Transform')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.plot(freq, amp)

# Fitting
if isFit:
    freq = freq[50:100]
    amp = amp[50:100]
    popt, pocv = curve_fit(gaussian, freq, amp, p0=[
                           1*10**3, 5.3*10**10, .2*10**11, 0])
    perr = np.sqrt(pocv)
    print(popt[0])
    print(popt[1])
    print(pocv)
    fit = gaussian(freq, *popt)
    plt.plot(freq, fit)
    plt.text(popt[1], popt[0]*1.25,
             "highest freq = {:e} +- {:e}".format(popt[1], popt[2]))


plt.xlim(-.1, plot_range[1])
plt.ylim(-200, 2*10**3)
plt.grid(False)
plt.savefig(outputFile)
