#ifndef HAMILTONDETECTOR_H
#define HAMILTONDETECTOR_H
#include <vector>
#include "RPeakDetector.h"

class HamiltonDetector : public RPeakDetector {
public:
    //refineRPeak = true: move each detection to the exact R-peak
    HamiltonDetector(int samplingRate, bool refineRPeak = true);
    std::vector<int> detect(const std::vector<double>& ecgSignal) override;

protected:
    //protected so the ablation variants can replace single stages
    int originalSamplingRate;
    int processingRate = 200;
    bool refineRPeak;

    //time windows in samples at 200 Hz
    int windowWidth = 17; //80 ms
    int ms95 = 19;
    int ms150 = 30;
    int ms200 = 40;
    int ms220 = 44;
    int ms360 = 72;
    int ms1000 = 200;
    int ms1500 = 300;

    double thresholdCoefficient = 0.475;
    int filterDelay = 59; //filter delay + 200 ms blanking
    int slopeBufferLength = 96;

    //signal processing
    std::vector<double> resample(const std::vector<double>& signal, int originalRate, int targetRate);
    static double sampleAt(const std::vector<double>& v, int i);
    virtual std::vector<double> lowPassFilter(const std::vector<double>& signal);
    virtual std::vector<double> highPassFilter(const std::vector<double>& signal);
    std::vector<double> derivative(const std::vector<double>& signal);
    std::vector<double> absoluteSignal(const std::vector<double>& signal);
    std::vector<double> movingWindowIntegration(const std::vector<double>& signal);

    //decision rules
    std::vector<int> findQRSComplexes(const std::vector<double>& integrated, const std::vector<double>& rawSlope);
    bool isBaselineShift(const std::vector<double>& rawSlope, int n, double& maxSlope);
    double findPeak(double value);
    double median(const std::vector<double>& values);
    virtual void updateThreshold();
    void acceptQRS(double peakHeight, int rrInterval);
    void restartLearning();
    virtual bool useSearchBack() const { return true; } //turned off in HamiltonAblation
    int findRPeak(const std::vector<double>& signal, int approximatePosition);

    //peak finder
    double peakMax = 0.0;
    int timeSinceMax = 0;
    double lastValue = 0.0;

    //200 ms blanking
    double heldPeak = 0.0;
    int blankCounter = 0;

    //last 8 values of each
    std::vector<double> qrsPeaks;
    std::vector<double> noisePeaks;
    std::vector<double> rrIntervals;

    double qrsMedian = 0.0;
    double noiseMedian = 0.0;
    double rrMedian = 0.0;
    double detectionThreshold = 0.0;

    int count = 0; //samples since last QRS
    int searchBackLimit = 0;

    //search-back candidate
    double searchBackPeak = 0.0;
    int searchBackLocation = 0;
    double searchBackSlope = 0.0;

    double lastQRSSlope = 0.0; //for the T-wave rule

    //learning phase and reset
    int learnedPeaks = 0;
    int secondCounter = 0;
    double secondMax = 0.0;
    std::vector<double> resetPeaks;

    std::vector<int> qrsPositions;
};

#endif