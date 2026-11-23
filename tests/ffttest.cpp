#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include "fft.h"
using namespace std;

//compares a transform with the textbook DFT formula and checks the inverse
void check(string name, vector<complex<double>> data, bool useFFT) {
    int n = data.size();
    vector<complex<double>> original = data;

    if (useFFT) {
        FFT::fft(data, false);
    } else {
        FFT::dft(data, false);
    }

    //slow textbook DFT for comparison
    double maxError = 0.0;
    for (int k = 0; k < n; k++) {
        complex<double> sum = 0;
        for (int t = 0; t < n; t++) {
            sum += original[t] * polar(1.0, -2 * M_PI * k * t / n);
        }
        maxError = max(maxError, abs(sum - data[k]));
    }

    //inverse should give back the original
    if (useFFT) {
        FFT::fft(data, true);
    } else {
        FFT::dft(data, true);
    }
    double backError = 0.0;
    for (int i = 0; i < n; i++) {
        backError = max(backError, abs(data[i] - original[i]));
    }

    cout << name << " (N = " << n << "): difference from textbook DFT = " << maxError
         << ", inverse error = " << backError << endl;
}

int main() {
    check("fft", {1, 2, 3, 4, 0, -1, -2, 5}, true); //power of 2
    check("dft", {1, 2, 3, 4, 0, -1, -2}, false); //N = 7, not a power of 2

    //longer test: N = 1000 (not a power of 2)
    vector<complex<double>> signal(1000);
    for (int i = 0; i < 1000; i++) {
        signal[i] = sin(0.3 * i) + 0.5 * cos(0.07 * i);
    }
    check("dft", signal, false);
}