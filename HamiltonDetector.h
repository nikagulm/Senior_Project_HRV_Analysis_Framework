#ifndef HAMILTONDETECTOR_H
#define HAMILTONDETECTOR_H
#include <vector>
#include "RPeakDetector.h"

class HamiltonDetector : public RPeakDetector {
    public:
        HamiltonDetector(int samplingRate);
        std::vector<int> detect(const std::vector<double>& ecgSignal) override;

    private:
        int samplingRate;
        std::vector<double> bandPassFilter(const std::vector<double>& signal);
        std::vector<double> derivative(const std::vector<double>& signal);
        std::vector<double> absoluteSignal(const std::vector<double>& signal);
        std::vector<double> movingAverage(const std::vector<double>& signal);
        std::vector<int> detectPeaks(const std::vector<double>& processedSignal, const std::vector<double>& originalSignal);

        int findRPeak(const std::vector<double>& signal, int candidatePosition);
    };

#endif