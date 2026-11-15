#ifndef FFTESTIMATOR_H
#define FFTESTIMATOR_H
#include "SpectrumEstimator.h"

class FFTEstimator : public SpectrumEstimator {
public:
    std::vector<double> power(const std::vector<double>& signal, double samplingRate, const std::vector<double>& frequencies) override;
};

#endif