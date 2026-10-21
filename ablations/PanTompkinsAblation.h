#ifndef PANTOMPKINSDETECTOR_H
#define PANTOMPKINSDETECTOR_H

#include <vector>
#include "RPeakDetector.h"

class PanTompkinsDetector : public RPeakDetector {
    public:
        PanTompkinsDetector(int samplingRate);
        std::vector<int> detect(const std::vector<double>& ecgSignal) override;

    private:
        int originalSamplingRate;
        int processingRate = 200;

        std::vector<double> resample(const std::vector<double>& signal, int originalRate, int targetRate);
        std::vector<double> bandPassFilter(const std::vector<double>& signal);
        std::vector<double> derivative(const std::vector<double>& signal);
        std::vector<double> squareSignal(const std::vector<double>& signal);
        std::vector<double> movingWindowIntegration(const std::vector<double>& signal);
        std::vector<int> detectPeaks(const std::vector<double>& integratedSignal, const std::vector<double>& originalSignal);
        int findRPeak(const std::vector<double>& signal,int candidatePosition);
    };

#endif