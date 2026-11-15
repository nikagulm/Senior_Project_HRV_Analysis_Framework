#include "FFTEstimator.h"
#include "fft.h"
#include <cmath>
#include <complex>
using namespace std;

vector<double> FFTEstimator::power(const vector<double>& signal, double samplingRate, const vector<double>& frequencies) {
    int m = paddedLength(signal.size());
    vector<complex<double>> data(m, 0.0);
    for (int i = 0; i < signal.size(); i++) {
        data[i] = signal[i];
    }

    //one FFT gives all m frequency components
    FFT::fft(data, false);
    vector<double> result;
    for (double f : frequencies) {
        int k = static_cast<int>(round(f * m / samplingRate));
        result.push_back(norm(data[k]));
    }
    return result;
}