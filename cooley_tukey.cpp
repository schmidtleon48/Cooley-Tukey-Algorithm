#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <fftw3.h>

using namespace std;

typedef complex<double> Complex;
typedef vector<Complex> CArray;

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
void bitReverseCopy(const CArray &a, CArray &A) {
    int n = a.size();
    int log2n = log2(n);
    for (int i = 0; i < n; i++) {
        int rev = reverseBits(i, log2n);
        A[rev] = a[i];
    }
}

// Iterative FFT function
CArray iterativeFFT(const CArray &a) {
    int n = a.size();
    CArray A(n);
    
    bitReverseCopy(a, A);  // Copy array with bit-reversed indexing
    
    for (int s = 1; s <= log2(n); s++) {
        int m = 1 << s;  // m = 2^s
        Complex wm = exp(Complex(0, -2.0 * M_PI / m));  // ωm = exp(−2πi/m)
        
        for (int k = 0; k < n; k += m) {
            Complex w = 1;
            for (int j = 0; j < m / 2; j++) {
                Complex t = w * A[k + j + m / 2];
                Complex u = A[k + j];
                A[k + j] = u + t;
                A[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }
    
    return A;
}

//Creates a vector with sine wave in real part and imaginary part zero
std::vector<fftw_complex> sinevector(int len_array, double start, double stop) {
    double step = (stop - start) / (len_array - 1);
    vector<fftw_complex> signal_sinus(len_array);

    for (int i = 0; i < len_array; i++) {
        double realPart = sin(static_cast<double>(i) * step);
        signal_sinus[i][0] = realPart;
        signal_sinus[i][1] = 0;
    }

    return signal_sinus;
}

// Main function
int main() {
    //Time series parameters
    const int len_array = 2048 * 4;
    const double start = 0.0;
    const double stop = 40.0 * 2.0 * M_PI;
    
    //Create time and frequency series
    std::vector<fftw_complex> sinus_time = sinevector(len_array, start, stop);
    std::vector<fftw_complex> sinus_frequency(len_array);

    // Create FFTW plan and execute it
    fftw_plan p = fftw_plan_dft_1d(len_array, sinus_time.data(), sinus_frequency.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    // Clean up FFTW resources
    fftw_destroy_plan(p);
    fftw_cleanup();

    return 0;
}