/*
Topic:          Cooley-Tukey-Algorithm
Description:    Compares FFTW and own implementation of Cooley-Tukey-Algorithm
Author:         Leon Schmidt
Date:           01.09.2024
Complile:       g++ cooley_tukey.cpp -o cooley_tukey -lfftw3 -lm
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include <chrono>

using namespace std;

// Function to reverse bits of an integer (used in bit-reverse-copy)
unsigned int reverseBits(unsigned int x, int log2n) {
    unsigned int n = 0;
    for (int i = 0; i < log2n; i++) {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }
    return n;
}

// Function to perform the bit-reverse copy
void bitReverseCopy(const vector<fftw_complex> &a, vector<fftw_complex> &A) {
    int n = a.size();
    int log2n = log2(n);
    for (int i = 0; i < n; i++) {
        int rev = reverseBits(i, log2n);
        A[rev][0] = a[i][0];
        A[rev][1] = a[i][1];
    }
}

// Iterative FFT function
vector<fftw_complex> iterativeFFT(const vector<fftw_complex> &a) {
    int n = a.size();
    vector<fftw_complex> A(n);
    
    bitReverseCopy(a, A);  // Copy array with bit-reversed indexing
    
    for (int s = 1; s <= log2(n); s++) {
        int m = 1 << s;  // m = 2^s
        fftw_complex wm;
        wm[0] = cos(-2.0 * M_PI / m);  // Real part of wm = cos(−2π/m)
        wm[1] = sin(-2.0 * M_PI / m);  // Imaginary part of wm = sin(−2π/m)
        
        for (int k = 0; k < n; k += m) {
            fftw_complex w = {1.0, 0.0};  // Initialize w to 1 (1 + 0i)
            for (int j = 0; j < m / 2; j++) {
                fftw_complex t = {
                    w[0] * A[k + j + m / 2][0] - w[1] * A[k + j + m / 2][1],  // Real part of t
                    w[0] * A[k + j + m / 2][1] + w[1] * A[k + j + m / 2][0]   // Imaginary part of t
                };
                fftw_complex u = {A[k + j][0], A[k + j][1]};  // Copy u = A[k + j]
                
                // Update A[k + j] and A[k + j + m / 2]
                A[k + j][0] = u[0] + t[0];  // Real part
                A[k + j][1] = u[1] + t[1];  // Imaginary part
                A[k + j + m / 2][0] = u[0] - t[0];  // Real part
                A[k + j + m / 2][1] = u[1] - t[1];  // Imaginary part

                // Update w: w *= wm
                double w_real = w[0] * wm[0] - w[1] * wm[1];
                double w_imag = w[0] * wm[1] + w[1] * wm[0];
                w[0] = w_real;
                w[1] = w_imag;
            }
        }
    }
    
    return A;
}

//Creates a vector with sine wave in real part and imaginary part zero, vector: len_array x 2
vector<fftw_complex> sinevector(int len_array, double start, double stop) {
    double step = (stop - start) / (len_array - 1);
    vector<fftw_complex> signal_sinus(len_array);

    for (int i = 0; i < len_array; i++) {
        double realPart = sin(static_cast<double>(i) * step);
        signal_sinus[i][0] = realPart;
        signal_sinus[i][1] = 0;
    }

    return signal_sinus;
}

// Function to calculate MSE for the real and imaginary parts
pair<double, double> calculateMSE(const vector<fftw_complex> &a, const vector<fftw_complex> &b) {
    if (a.size() != b.size()) {
        cerr << "Vectors must be of the same size." << endl;
        return {-1.0, -1.0};  // Return an error indicator
    }

    double mse_real = 0.0;
    double mse_imag = 0.0;
    int n = a.size();

    for (int i = 0; i < n; ++i) {
        double diff_real = a[i][0] - b[i][0];
        double diff_imag = a[i][1] - b[i][1];

        mse_real += diff_real * diff_real;
        mse_imag += diff_imag * diff_imag;
    }

    mse_real /= n;
    mse_imag /= n;

    return {mse_real, mse_imag};  
}

// Main function
int main() {
    // Time series parameters
    const int len_array = static_cast<int>(pow(2, 20));
    const double start = 0.0;
    const double stop = 40.0 * 2.0 * M_PI;
    
    // Create time and frequency series
    vector<fftw_complex> sinus_time = sinevector(len_array, start, stop);
    vector<fftw_complex> sinus_frequency_fftw(len_array);

    // Measure performance of FFTW
    auto start_fftw = chrono::high_resolution_clock::now();

    // Create FFTW plan and execute it
    fftw_plan p = fftw_plan_dft_1d(len_array, sinus_time.data(), sinus_frequency_fftw.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Measure end time of FFTW
    auto end_fftw = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_fftw = end_fftw - start_fftw;

    // Clean up FFTW resources
    fftw_destroy_plan(p);
    fftw_cleanup();

    // Measure performance of own iterative FFT implementation
    auto start_iterative = chrono::high_resolution_clock::now();
    
    vector<fftw_complex> sinus_frequency_iterative = iterativeFFT(sinus_time);

    auto end_iterative = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_iterative = end_iterative - start_iterative;

    // Calculate MSE
    pair<double, double> mse = calculateMSE(sinus_frequency_fftw, sinus_frequency_iterative);

    // Output the length of the array
    cout << "Length array in log2: " << log2(len_array) << "\n";

    // Output the MSE results
    cout << "MSE for real part: " << mse.first << endl;
    cout << "MSE for imaginary part: " << mse.second << "\n";

    // Output the timing results
    cout << "Time taken by FFTW: " << elapsed_fftw.count() << " seconds" << endl;
    cout << "Time taken by iterative FFT: " << elapsed_iterative.count() << " seconds" << endl;

    return 0;
}