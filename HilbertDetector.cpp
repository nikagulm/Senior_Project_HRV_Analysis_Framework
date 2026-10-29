#include "HilbertDetector.h"
#include "fft.h"
#include <algorithm>
#include <cmath>
using namespace std;

HilbertDetector::HilbertDetector(int samplingRate) {
    originalSamplingRate = samplingRate;
}


//main function to return the r-peaks
vector<int> HilbertDetector::detect(const vector<double>& ecgSignal) {
    if (ecgSignal.empty()) {
        return {};
    }

    //I.resample to 360 Hz cuz of mitbih
    vector<double> signal = ecgSignal;
    if (originalSamplingRate != processingRate) {
        signal = resample(ecgSignal, originalSamplingRate, processingRate);
    }

    //II.bandpass 8-20Hz then derivative
    vector<double> filtered = bandPassFilter(signal);
    vector<double> derivativeSignal = derivative(filtered);

    //III. hilbert transform of the derivative
    vector<double> hilbert = hilbertTransform(derivativeSignal);

    //IV: threshold for each segment + find peaks above it
    vector<double> thresholds = segmentThresholds(hilbert);
    vector<int> peaks = findPeaks(hilbert, thresholds);

    //V: find the actual r peak in the ecg and go back to the original rate
    vector<int> rPeaks;
    for (int peak : peaks) {
        int rPeak = findRPeak(signal, peak);
        int originalPosition = static_cast<int>(round(static_cast<double>(rPeak) * originalSamplingRate / processingRate));
        if (rPeaks.empty() || originalPosition > rPeaks.back()) {
            rPeaks.push_back(originalPosition);
        }
    }
    return rPeaks;
}



//linear interpolation > same as in the other detectors
vector<double> HilbertDetector::resample(const vector<double>& signal, int originalRate, int targetRate) {
    int newSize = static_cast<int>(static_cast<double>(signal.size()) * targetRate / originalRate);
    vector<double> result(newSize);
    for (int i = 0; i < newSize; i++) {
        double originalPosition = static_cast<double>(i) * originalRate / targetRate;
        int left = static_cast<int>(originalPosition);
        int right = min(left + 1, static_cast<int>(signal.size()) - 1);
        double fraction = originalPosition - left;
        result[i] = signal[left] * (1.0 - fraction) + signal[right] * fraction;
    }
    return result;
}



//bandpass filter 8-20Hz with a kaiser window
//not speicified in original alg so picked 60dB and 4Hz transition -> about 327 coefficients
//also shifted the output so theres no delay
vector<double> HilbertDetector::bandPassFilter(const vector<double>& signal) {
    double fs = processingRate;
    double lowCut = 8.0;
    double highCut = 20.0;

    //kaiser stuff (formulas for beta and order)
    double attenuation = 60.0;
    double transition = 4.0;
    double beta = 0.1102 * (attenuation - 8.7);
    int order = static_cast<int>(ceil((attenuation - 8.0) / (2.285 * 2.0 * M_PI * transition / fs)));
    if (order % 2 == 1) {
        order++;
    }
    int length = order + 1;
    int middle = order / 2;

    //bessel function I0, the kaiser window needs it
    auto besselI0 = [](double x) {
        double sum = 1.0;
        double term = 1.0;
        for (int k = 1; k < 50; k++) {
            term *= (x / (2.0 * k)) * (x / (2.0 * k));
            sum += term;
        }
        return sum;
    };

    //ideal bandpass times the window = filter coefficients
    vector<double> h(length);
    for (int n = 0; n < length; n++) {
        int k = n - middle;
        double ideal;
        if (k == 0) {
            ideal = 2.0 * (highCut - lowCut) / fs;
        } else {
            ideal = (sin(2.0 * M_PI * highCut * k / fs) - sin(2.0 * M_PI * lowCut * k / fs)) / (M_PI * k);
        }

        double r = static_cast<double>(k) / middle;
        double window = besselI0(beta * sqrt(1.0 - r * r)) / besselI0(beta);
        h[n] = ideal * window;
    }

    //convolution (shifted by middle so it lines up with the input)
    int size = signal.size();
    vector<double> output(size, 0.0);
    for (int i = 0; i < size; i++) {
        double sum = 0.0;
        for (int n = 0; n < length; n++) {
            int j = i + middle - n;
            if (j >= 0 && j < size) {
                sum += h[n] * signal[j];
            }
        }
        output[i] = sum;
    }
    return output;
}



//derivative y(n) = 0.5*(x(n+1) - x(n-1))
vector<double> HilbertDetector::derivative(const vector<double>& x) {
    int size = x.size();
    vector<double> y(size, 0.0);

    for (int n = 1; n < size - 1; n++) {
        y[n] = 0.5 * (x[n + 1] - x[n - 1]);
    }
    return y;
}


//hilbert transform on the whole signal (did it per segment first but it gave false peaks at the edges)
//dft -> positive freqs * -j, negative freqs * +j -> inverse dft -> take the real part
vector<double> HilbertDetector::hilbertTransform(const vector<double>& signal) {
    int size = signal.size();

    vector<complex<double>> data(size);
    for (int i = 0; i < size; i++) {
        data[i] = signal[i];
    }

    FFT::dft(data, false);

    //the 90 degree shift n 0Hz and the middle bin just become 0
    data[0] = 0.0;
    for (int k = 1; k < (size + 1) / 2; k++) {
        data[k] *= complex<double>(0.0, -1.0);
        data[size - k] *= complex<double>(0.0, 1.0);
    }
    if (size % 2 == 0) {
        data[size / 2] = 0.0;
    }

    FFT::dft(data, true);
    vector<double> result(size);
    for (int i = 0; i < size; i++) {
        result[i] = data[i].real();
    }

    return result;
}


//threshold for every 1024 samples (eq 4 in arzeno)
//rms is 0.18*max -> 0.39*max, or 0.39*previous max if max jumped more than 2x
//else its 1.6*rms
vector<double> HilbertDetector::segmentThresholds(const vector<double>& hilbert) {
    int size = hilbert.size();
    vector<double> thresholds(size, 0.0);
    double previousMax = -1.0; //-1 = first segment cuz no previous yet

    for (int start = 0; start < size; start += segmentLength) {
        int end = min(start + segmentLength, size);
        double maxValue = 0.0;
        double sumSquares = 0.0;
        for (int i = start; i < end; i++) {
            maxValue = max(maxValue, hilbert[i]);
            sumSquares += hilbert[i] * hilbert[i];
        }
        double rms = sqrt(sumSquares / (end - start));

        double threshold;
        if (rms > 0.18 * maxValue) {
            //big jump is probably noise, so use the previous max
            //(skipped if previous max is 0, otherwise threshold would be 0)
            if (previousMax > 0 && maxValue > 2.0 * previousMax) {
                threshold = 0.39 * previousMax;
            } else {
                threshold = 0.39 * maxValue;
            }
        } else {
            threshold = 1.6 * rms;
        }

        for (int i = start; i < end; i++) {
            thresholds[i] = threshold;
        }
        previousMax = maxValue;
    }
    return thresholds;
}



//when it goes above the threshold take the biggest value in the next 200ms
//then skip 200ms after that peak
vector<int> HilbertDetector::findPeaks(const vector<double>& hilbert, const vector<double>& thresholds) {
    vector<int> peaks;
    int size = hilbert.size();
    int n = 0;
    while (n < size) {
        if (hilbert[n] > thresholds[n]) {
            int peak = n;
            int end = min(n + refractory, size);
            for (int i = n; i < end; i++) {
                if (hilbert[i] > hilbert[peak]) {
                    peak = i;
                }
            }

            peaks.push_back(peak);
            n = peak + refractory; //refractory period
        } else {
            n++;
        }
    }
    return peaks;
}



//r peak = max of the ecg within +-10 samples of the hilbert peak
int HilbertDetector::findRPeak(const vector<double>& signal, int position) {
    int start = max(0, position - searchRadius);
    int end = min(static_cast<int>(signal.size()) - 1, position + searchRadius);
    int rPeak = start;
    for (int i = start + 1; i <= end; i++) {
        if (signal[i] > signal[rPeak]) {
            rPeak = i;
        }
    }
    return rPeak;
}