import numpy as np
import matplotlib.pyplot as plt

inputFile = "../dat/theoretical/ref/lambda/g_main.dat"

# Load data from file
data = np.loadtxt(inputFile, delimiter=' ')
def eq(x, y): return (x-y)/(x+y)

# Extract time and amplitude data
time = list(data[:,0])
amplitude = list(map(eq, data[:,1], data[:,2]))

# Compute FFT
fft_result = np.fft.fft(amplitude)
# Compute corresponding frequencies
frequency = np.fft.fftfreq(len(time), np.mean(np.diff(time)))

# Plot magnitude of FFT
plt.figure(figsize=(8, 5))
plt.title('Fourier Transform')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.plot(np.abs(frequency), np.abs(fft_result))
plt.grid(False)
plt.show()
