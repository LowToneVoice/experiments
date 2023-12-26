import numpy as np
import matplotlib.pyplot as plt

inputFile = "../dat/theoretical/ref/lambda/g_main.dat"

# Load data from file
data = np.loadtxt(inputFile, delimiter=' ')
def eq(x, y): return (x-y)/(x+y)

wavelength=list(data[:,0])
data = list(map(eq, data[:, 1], data[:,2]))

print(wavelength)

