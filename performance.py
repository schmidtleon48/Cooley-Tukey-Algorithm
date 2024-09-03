'''
Topic:          Cooley-Tukey-Algorithm Performance and Validation
Description:    Compares FFTW and own implementation of Cooley-Tukey-Algorithm, Performance and Validation
Author:         Leon Schmidt
Date:           01.09.2024
'''

import matplotlib.pyplot as plt
import numpy as np

len_array = np.linspace(2**10, 2**24, 1000)

# The values are gained by performing the cooley_tukey algorithm with different array lengths
len_array_test = np.array([2**10, 2**12, 2**16, 2**20, 2**24])
time_fftw = np.array([0.00064, 0.00066, 0.0039, 0.043, 1.09])
time_iterative = np.array([0.00022, 0.001, 0.016, 0.383, 7.37])

faktor = (time_iterative[0]) / (len_array[0] * np.log2(len_array[0]))
plt.scatter(len_array_test, time_fftw, label = "FFT in the West")
plt.scatter(len_array_test, time_iterative, label = "Iterative FFT")
plt.plot(len_array, faktor * len_array * np.log2(len_array), linestyle='--', label = "O(n * log2(n))")
plt.xscale('log', base=2)
plt.yscale('log', base=10)
plt.title("Performance Analysis")
plt.xlabel("Length Time Series")
plt.ylabel("Executing Time [s]")
plt.legend()
plt.plot()
plt.savefig("images/performance.png")
plt.show()
plt.clf()


mse_real = np.array([4e-27, 1.25e-25, 1.25e-21, 6.5e-18, 2.7e-14])
mse_imag = np.array([5.3e-26, 2.46e-24, 1.38e-20, 5.4e-17, 3.16e-13])

plt.plot(len_array_test, mse_real, linestyle = '--', marker = "x", label = "MSE Real part")
plt.plot(len_array_test, mse_imag, linestyle = '--', marker = "x", label = "MSE Imag part")
plt.xscale('log', base=2)
plt.yscale('log', base=10)
plt.title("Mean Squared Error")
plt.xlabel("Length Time Series")
plt.ylabel("MSE")
plt.legend()
plt.plot()
plt.savefig("images/mse.png")
plt.show()