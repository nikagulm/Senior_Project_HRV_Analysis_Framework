#include "HamiltonAblation.h"
using namespace std;

HamiltonAblation::HamiltonAblation(int samplingRate, AblationVariant variant)
    : HamiltonDetector(samplingRate), variant(variant) {
    if (variant == AblationVariant::MovingAverage) {
        //take out the low-pass (4) and high-pass (12.5) delays, add the moving average delay
        filterDelay = filterDelay - 4 - 12 + (averageLength - 1) / 2;
    }
}



vector<int> HamiltonAblation::detect(const vector<double>& ecgSignal) {
    thresholdFrozen = false;
    return HamiltonDetector::detect(ecgSignal);
}


//moving average variant: the low-pass becomes a simple moving average
vector<double> HamiltonAblation::lowPassFilter(const vector<double>& x) {
    if (variant != AblationVariant::MovingAverage) {
        return HamiltonDetector::lowPassFilter(x);
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
vector<double> HamiltonAblation::highPassFilter(const vector<double>& x) {
    if (variant != AblationVariant::MovingAverage) {
        return HamiltonDetector::highPassFilter(x);
    }
    return x;
}

//fixed threshold variant: the first threshold after learning is kept for the whole record
void HamiltonAblation::updateThreshold() {
    if (thresholdFrozen) {
        return;
    }
    HamiltonDetector::updateThreshold();
    if (variant == AblationVariant::FixedThreshold && learnedPeaks == 8) {
        thresholdFrozen = true;
    }
}



//no search-back variant
bool HamiltonAblation::useSearchBack() const {
    return variant != AblationVariant::NoSearchBack;
}