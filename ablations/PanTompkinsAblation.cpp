#include "PanTompkinsAblation.h"
using namespace std;

PanTompkinsAblation::PanTompkinsAblation(int samplingRate, AblationVariant variant)
    : PanTompkinsDetector(samplingRate), variant(variant) {
    if (variant == AblationVariant::MovingAverage) {
        //delay of a centred moving average is (length - 1) / 2, and there is no high-pass
        lowPassDelay = (averageLength - 1) / 2;
        highPassDelay = 0;
    }
}



vector<int> PanTompkinsAblation::detect(const vector<double>& ecgSignal) {
    thresholdsFrozen = false;
    return PanTompkinsDetector::detect(ecgSignal);
}



//moving average variant: the low-pass becomes a simple moving average
vector<double> PanTompkinsAblation::lowPassFilter(const vector<double>& x) {
    if (variant != AblationVariant::MovingAverage) {
        return PanTompkinsDetector::lowPassFilter(x);
    }

    int size = x.size();
    vector<double> y(size, 0.0);
    double sum = 0.0;

    for (int n = 0; n < size; n++) {
        sum += x[n];
        if (n >= averageLength) {
            sum -= x[n - averageLength];
        }
        y[n] = sum / averageLength;
    }

    return y;
}



//moving average variant: no high-pass, so baseline wander stays in the signal
vector<double> PanTompkinsAblation::highPassFilter(const vector<double>& x) {
    if (variant != AblationVariant::MovingAverage) {
        return PanTompkinsDetector::highPassFilter(x);
    }
    return x;
}



//fixed threshold variant: learn as usual, then freeze the thresholds
int PanTompkinsAblation::learnThresholds(const vector<double>& integrated, const vector<double>& bandPassed) {
    int end = PanTompkinsDetector::learnThresholds(integrated, bandPassed);
    if (variant == AblationVariant::FixedThreshold) {
        thresholdsFrozen = true;
    }
    return end;
}



void PanTompkinsAblation::updateThresholds() {
    if (thresholdsFrozen) {
        return;
    }
    PanTompkinsDetector::updateThresholds();
}



//no search-back variant: missed beats are never recovered
bool PanTompkinsAblation::searchBack() {
    if (variant == AblationVariant::NoSearchBack) {
        return false;
    }
    return PanTompkinsDetector::searchBack();
}