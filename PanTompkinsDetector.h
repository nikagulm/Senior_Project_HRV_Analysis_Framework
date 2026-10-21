#ifndef PANTOMPKINSDETECTOR_H
#define PANTOMPKINSDETECTOR_H
#include <vector>
#include "RPeakDetector.h"


// one peak found in the integrated signal
struct PeakCandidate {
    int position; // position in the integrated signal
    double peakI; // height of the peak in the integrated signal (PEAKI)
    double peakF; // height of the matching peak in the band-passed signal (PEAKF)
    int positionF; // where PEAKF is in the band-passed signal
    double slope; // steepest slope of the band-passed signal (for the T-wave check)
};

class PanTompkinsDetector : public RPeakDetector {
public:
    PanTompkinsDetector(int samplingRate);

    std::vector<int> detect(const std::vector<double>& ecgSignal) override;

private:
    int originalSamplingRate;
    int processingRate = 200;

    // delays of the filters in samples (at 200 Hz)
    int lowPassDelay = 5;
    int highPassDelay = 16;

    //signal processing stages
    std::vector<double> resample(const std::vector<double>& signal, int originalRate, int targetRate);
    static double sampleAt(const std::vector<double>& v, int i);
    std::vector<double> lowPassFilter(const std::vector<double>& signal);
    std::vector<double> highPassFilter(const std::vector<double>& signal);
    std::vector<double> derivative(const std::vector<double>& signal);
    std::vector<double> squareSignal(const std::vector<double>& signal);
    std::vector<double> movingWindowIntegration(const std::vector<double>& signal);

    //decision stage 
    std::vector<int> findQRSComplexes(const std::vector<double>& integrated, const std::vector<double>& bandPassed);

    PeakCandidate makeCandidate(int position, const std::vector<double>& integrated, const std::vector<double>& bandPassed);
    int learnThresholds(const std::vector<double>& integrated, const std::vector<double>& bandPassed);
    void updateThresholds();
    void updateRRAverages(int rrInterval);
    bool isTWave(const PeakCandidate& candidate);
    void acceptQRS(const PeakCandidate& candidate);
    void markAsNoise(const PeakCandidate& candidate);
    bool searchBack();
    int findRPeak(const std::vector<double>& signal, int approximatePosition);

    //thresholds for the integrated signal
    double signalPeakI = 0.0;// SPKI
    double noisePeakI = 0.0; // NPKI
    double thresholdI1 = 0.0; // THRESHOLD I1
    double thresholdI2 = 0.0; // THRESHOLD I2 (used in search-back)

    //thresholds for the band-passed signal
    double signalPeakF = 0.0; // SPKF
    double noisePeakF = 0.0;  // NPKF
    double thresholdF1 = 0.0; // THRESHOLD F1
    double thresholdF2 = 0.0;  // THRESHOLD F2 (used in search-back)

    //RR-interval information
    std::vector<int> recentRR; // last 8 RR intervals
    std::vector<int> regularRR; // last 8 RR intervals inside the 92%-116% limits
    double rrAverage1 = 0.0;  // RR AVERAGE1
    double rrAverage2 = 0.0;   // RR AVERAGE2
    double rrMissedLimit = 0.0;   // RR MISSED LIMIT
    bool irregularRhythm = false;

    //detection state 
    int lastQRSPosition = -1;
    double lastQRSSlope = 0.0;
    std::vector<PeakCandidate> noisePeaks;   // noise peaks since the last QRS
    std::vector<int> qrsPositions; //detected QRS positions (band-passed signal)
};

#endif