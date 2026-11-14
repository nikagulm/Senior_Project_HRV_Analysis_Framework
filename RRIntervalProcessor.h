#ifndef RRINTERVALPROCESSOR_H
#define RRINTERVALPROCESSOR_H
#include <vector>

//turns R-peak positions into RR intervals and the cleaned NN sequence
class RRIntervalProcessor {
public:
    RRIntervalProcessor(int samplingRate);
    std::vector<double> rrIntervals(const std::vector<int>& rPeaks);
    std::vector<double> nnIntervals(const std::vector<double>& rr);

private:
    int samplingRate;
    double minInterval = 0.3; //200 bpm
    double maxInterval = 2.0; //30 bpm
    double maxChange = 0.2; //20%
};

#endif