#ifndef FREQUENCYDOMAINHRV_H
#define FREQUENCYDOMAINHRV_H
#include <vector>
#include "SpectrumEstimator.h"

struct FrequencyDomainResult {
    double vlf;
    double lf;
    double hf;
    double lfhf;
    int segments; 
};

class FrequencyDomainHRV {
public:
    FrequencyDomainHRV(SpectrumEstimator* estimator);
    FrequencyDomainResult calculate(const std::vector<double>& nnTimes, const std::vector<double>& nn);

private:
    SpectrumEstimator* estimator;
    double resampleRate = 4.0; //Hz
    int segmentLength = 1200;  //5 minutes at 4 Hz
    
    //band limits in hz
    double vlfLow = 0.0033;
    double lfLow = 0.04;
    double hfLow = 0.15;
    double hfHigh = 0.40;
};

#endif