#ifndef HILBERTDETECTOR_H
#define HILBERTDETECTOR_H
#include <vector>
#include "RPeakDetector.h"



class HilbertDetector : public RPeakDetector {
public:
    HilbertDetector(int samplingRate);
    std::vector<int> detect(const std::vector<double>& ecgSignal) override;

private:
    int originalSamplingRate;
    int processingRate = 360; //the algorithm was designed for MIT-BIH data at 360 Hz
    int segmentLength = 1024; //samples per threshold window
    int refractory = 72; //200 ms at 360 Hz
    int searchRadius = 10; //R-peak search +-10 samples around the Hilbert peak

    //signal processing stages
    std::vector<double> resample(const std::vector<double>& signal, int originalRate, int targetRate);
    std::vector<double> bandPassFilter(const std::vector<double>& signal);
    std::vector<double> derivative(const std::vector<double>& signal);
    std::vector<double> hilbertTransform(const std::vector<double>& signal);

    //decision stage
    std::vector<double> segmentThresholds(const std::vector<double>& hilbert);
    std::vector<int> findPeaks(const std::vector<double>& hilbert, const std::vector<double>& thresholds);
    int findRPeak(const std::vector<double>& signal, int position);
};

#endif