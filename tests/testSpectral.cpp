//Fourier transforms, Goertzel and frequency-domain HRV
#include "TestUtils.h"
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include "fft.h"
#include "FFTEstimator.h"
#include "GoertzelEstimator.h"
#include "FrequencyDomainHRV.h"
using namespace std;

//DFT formula for reference
vector<complex<double>> slowDFT(const vector<complex<double>>& x) {
    int n = x.size();
    vector<complex<double>> X(n);
    for (int k = 0; k < n; k++) {
        complex<double> sum = 0;
        for (int t = 0; t < n; t++) {
            sum += x[t] * polar(1.0, -2.0 * M_PI * k * t / n);
        }
        X[k] = sum;
    }
    return X;
}

//largest difference between two spectra relative to the largest value
double relativeError(const vector<complex<double>>& a, const vector<complex<double>>& b) {
    double err = 0, top = 0;
    for (int i = 0; i < a.size(); i++) {
        err = max(err, abs(a[i] - b[i]));
        top = max(top, abs(b[i]));
    }
    return err / top;
}

//spectral analysis
void testTransforms() {
    printf("Fourier transforms");
    mt19937 rng(1);
    normal_distribution<double> noise(0.0, 1.0);
    for (int n : {8, 64, 1024, 2048, 7, 100, 1000, 1201}) {
        vector<complex<double>> x(n);
        for (auto& v : x) { 
            v = {noise(rng), noise(rng)};
        }
        bool powerOfTwo = (n & (n - 1)) == 0;
        vector<complex<double>> X = x;
        if (powerOfTwo) {
            FFT::fft(X, false);
        } else {
            FFT::dft(X, false);
        }
        double e = relativeError(X, slowDFT(x));
        check(string(powerOfTwo ? "radix-2 FFT" : "Bluestein DFT") + " vs textbook DFT, N = " + to_string(n), e < 1e-9, fmt("relative error %.1e", e));

        //inverse transform must give the input back
        if (powerOfTwo) {
            FFT::fft(X, true);
        } else {
            FFT::dft(X, true);
        }
        double r = relativeError(X, x);
        check(" inverse transform returns the input, N = " + to_string(n), r < 1e-9, fmt("relative error %.1e", r));
    }
}

void testGoertzel() {
    printf("Goertzel vs FFT ");
    mt19937 rng(2);
    normal_distribution<double> noise(0.0, 30.0);
    vector<double> segment(1200);
    for (double& v : segment) {
        v = noise(rng);
    }
    vector<double> frequencies;
    for (int k = 1; k * 4.0 / 2048 <= 0.40; k++) {
        frequencies.push_back(k * 4.0 / 2048);
    }
    FFTEstimator fft;
    GoertzelEstimator goertzel;
    vector<double> a = fft.power(segment, 4.0, frequencies);
    vector<double> b = goertzel.power(segment, 4.0, frequencies);
    double err = 0, top = 0;
    for (int i = 0; i < a.size(); i++) {
        err = max(err, abs(a[i] - b[i]));
        top = max(top, abs(a[i]));
    }
    check("Goertzel equals FFT on " + to_string(frequencies.size()) + " frequencies (1200 samples)", err / top < 1e-9, fmt("relative difference %.1e", err / top));
}

//NN series made of known waves sampled every 0.25s for 25 minutes
FrequencyDomainResult spectrumOf(double vlfAmp, double lfAmp, double hfAmp) {
    vector<double> times, nn;
    for (int i = 1; i <= 6000; i++) {
        double t = i * 0.25;
        times.push_back(t);
        double ms = 800 + vlfAmp * sin(2 * M_PI * 0.03 * t) + lfAmp * sin(2 * M_PI * 0.10 * t) + hfAmp * sin(2 * M_PI * 0.25 * t);
        nn.push_back(ms / 1000.0);
    }
    FFTEstimator fft;
    FrequencyDomainHRV fd(&fft);
    return fd.calculate(times, nn);
}

void testBandPowers() {
    printf("Frequency-domain HRV ");
    FrequencyDomainResult r = spectrumOf(0, 30, 20);
    check("LF power of a 30 ms wave at 0.10 Hz (expected 450 ms^2)", abs(r.lf - 450) < 0.05 * 450, fmt("LF = %.1f ms^2", r.lf));
    check("HF power of a 20 ms wave at 0.25 Hz (expected 200 ms^2)", abs(r.hf - 200) < 0.05 * 200, fmt("HF = %.1f ms^2", r.hf));
    check("LF/HF ratio (expected 2.25)", abs(r.lfhf - 2.25) < 0.1, fmt("LF/HF = %.3f", r.lfhf));
    check("number of five-minute segments in 25 minutes (expected 5)", r.segments == 5, fmt("%.0f segments", r.segments));

    //a VLF wave must stay in VLF the Hann window keeps its power out of LF
    FrequencyDomainResult v = spectrumOf(30, 0, 0);
    check("0.03 Hz wave: power leaking into LF (Hann window)", v.lf < 1.0, fmt("VLF = %.1f ms^2, LF = %.2f ms^2", v.vlf, v.lf));
}