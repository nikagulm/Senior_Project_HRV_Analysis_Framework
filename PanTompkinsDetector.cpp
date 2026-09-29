#include "PanTompkinsDetector.h"
using namespace std;

vector<int> PanTompkinsDetector::detect(const vector<double>& ecgSignal) {
    return {};
}

//applies the low-pass stage of the Pan-Tompkins preprocessing
vector<double> PanTompkinsDetector::bandPassFilter(const vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }
    vector<double> filteredValues(signal.size(), 0);
    //low-pass filter
    for (int i = 12; i < signal.size(); i++) {
        filteredValues[i] =
                2 * filteredValues[i - 1]
                - filteredValues[i - 2]
                + signal[i]
                - 2 * signal[i - 6]
                + signal[i - 12];
    }
    return filteredValues;
}

//tells where signal changes quickly
vector<double> PanTompkinsDetector::derivative(const vector<double>& signal) {
    vector <double> derivativeValues;
    for (int i = 1; i < signal.size(); i++) {
        derivativeValues.push_back(signal[i] - signal[i - 1]);
}
    return derivativeValues;
}

//squares every value in the signal so large changes become more noticeable
vector<double> PanTompkinsDetector::squareSignal(const vector<double>& signal) {
    vector<double> squareValues;
    for (double value : signal) {
        squareValues.push_back(value*value);
    }
    return squareValues;
}

vector<double> PanTompkinsDetector::movingWindowIntegration(const vector<double>& signal) {
    return {};
}

vector<int> PanTompkinsDetector::detectPeaks(const vector<double>& signal) {
    return {};
}