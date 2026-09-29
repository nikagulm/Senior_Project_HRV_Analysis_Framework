#include "PanTompkinsDetector.h"
using namespace std;

vector<int> PanTompkinsDetector::detect(const vector<double>& ecgSignal) {
    return {};
}

//applies the low and high-pass stage of the Pan-Tompkins preprocessing for very fast and slow noise
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
    
    vector<double> highPassValues(signal.size(), 0);
    //high-pass filter
    for (int i = 32; i < signal.size(); i++) {
        highPassValues[i] =
                highPassValues[i - 1]
                - filteredValues[i] / 32.0
                + filteredValues[i - 16]
                - filteredValues[i - 17]
                + filteredValues[i - 32] / 32.0;
    }
    return highPassValues;
}


//tells where signal changes quickly
vector<double> PanTompkinsDetector::derivative(const vector<double>& signal) {
    if (signal.size() < 5) {
        return {};
    }
    vector<double> derivativeValues(signal.size(), 0);
    for (int i = 2; i < signal.size() - 2; i++) {
        derivativeValues[i] =
                (-signal[i - 2]
                 - 2 * signal[i - 1]
                 + 2 * signal[i + 1]
                 + signal[i + 2]) / 8.0;
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

//averages the squared signal over a short window to emphasize QRS complexes
vector<double> PanTompkinsDetector::movingWindowIntegration(const vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }
    vector<double> integratedValues(signal.size(), 0);
    //150 ms window for a 200 Hz signal
    int windowSize = 30;

    for (int i = windowSize - 1; i < signal.size(); i++) {
        double sum = 0;
        //add all values inside the current window
        for (int j = i - windowSize + 1; j <= i; j++) {
            sum += signal[j];
        }
        //store the average of the window
        integratedValues[i] = sum / windowSize;
    }
    return integratedValues;
}

//detects QRS peak positions from the processed ECG signal
vector<int> PanTompkinsDetector::detectPeaks(const vector<double>& signal) {
    vector<int> peaks;
    if (signal.empty()) {
        return peaks;
    }
    
    double noisePeakEstimate = 0.0; 
    double threshold = 0.0;
    double signalPeakEstimate = 0.0;

    int lastPeak = -1; //last heartbeat received to check for fake ones
    vector <int> rrIntervals; //distance between peaks
    int refractoryPeriod = 40;

    for (int i = 1; i<signal.size() - 1; i++) {
        if (signal[i] > signal[i - 1] && signal[i] > signal[i+1]){
            double peakValue = signal[i];
            if (peakValue > threshold){
                if (lastPeak == -1 || i - lastPeak > refractoryPeriod){
                    signalPeakEstimate = 0.125 * peakValue + 0.875 * signalPeakEstimate;
                
                    peaks.push_back(i); //i is the detected peak
                    lastPeak = i;
            } else {
                noisePeakEstimate = 0.125 * peakValue + 0.875 * noisePeakEstimate;
            }
            threshold = noisePeakEstimate + 0.25 * (signalPeakEstimate - noisePeakEstimate);
        }
    }
    return peaks;
}



