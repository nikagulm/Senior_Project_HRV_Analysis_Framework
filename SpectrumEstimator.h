#ifndef SPECTRUMESTIMATOR_H
#define SPECTRUMESTIMATOR_H
#include <vector>

//common interface for FFT and Goertzel
class SpectrumEstimator {
public:
    virtual ~SpectrumEstimator() = default;
    virtual std::vector<double> power(const std::vector<double>& signal, double samplingRate, const std::vector<double>& frequencies) = 0;

protected:
    static int paddedLength(int n) {
        int m = 1;
        while (m < n) {
            m *= 2;
        }
        return m;
    }
};

#endif