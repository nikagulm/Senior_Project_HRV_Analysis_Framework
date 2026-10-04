#include "PanTompkinsDetector.h"
using namespace std;


PanTompkinsDetector::PanTompkinsDetector(int samplingRate) {
    this->samplingRate = samplingRate;
}

vector<int> PanTompkinsDetector::detect(const vector<double>& ecgSignal) {

    vector<double> filteredSignal = bandPassFilter(ecgSignal);
    vector<double> derivativeSignal = derivative(filteredSignal);
    vector<double> squaredSignal = squareSignal(derivativeSignal);
    vector<double> integratedSignal =
        movingWindowIntegration(squaredSignal);
    vector<int> peaks = detectPeaks(integratedSignal);

    return peaks;
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

    for (int i = windowSize -1; i < signal.size();i++) {
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

//cetects QRS peak positions from the processed ECG signal
vector<int> PanTompkinsDetector::detectPeaks(const vector<double>& signal) {
    vector<int> peaks;
    if (signal.empty()) {
        return peaks;
    }

    //for choosing reasonable starting threshold 
    //fix later 
    int initializationSamples = 2 * samplingRate; 
    if (signal.size() < initializationSamples) {
        initializationSamples = signal.size();
    }
    double maxValue = signal[0];
    double sumValue = 0.0;

    for (int i = 0; i < initializationSamples; i++) {
        sumValue += signal[i];
        if (signal[i] > maxValue) {
         maxValue = signal[i];
        }
    }

    double averageValue = sumValue / initializationSamples;
    // Set the starting estimates before peak detection begins
    //fix THE BELOW LATER
    double signalPeakEstimate = 0.25 * maxValue;
    double noisePeakEstimate = 0.5 * averageValue;
    double threshold = noisePeakEstimate + 0.25 * (signalPeakEstimate - noisePeakEstimate);
    double searchBackThreshold = 0.5 * threshold;

    //threshold values -> TO BE UPDATED/REMOVED 
    //double noisePeakEstimate = 0.0;
    //double threshold = 0.0;
    //double signalPeakEstimate = 0.0;
    //double searchBackThreshold = 0.0;
   
    //these in case for preventing counting one heartbeat twice
    int lastPeak = -1;
    int refractoryPeriod = 0.2 * samplingRate; //period for minimum distance bw heartbeats
    int tWaveWindow = 0.36 * samplingRate; //time window where new peak may be t wave

    vector<int> rrIntervals; //distance between peaks
    double rrAverage = 0.0;
    double missedBeatLimit = 0.0;

    //for storing the best weak peak in case a heartbeat was missed [search-back process]
    int searchBackPeak = -1;
    double searchBackValue = 0.0;

    for (int i = 1; i < signal.size() - 1; i++) {
    //check if current point is a local maximum
        if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1]) {
            double peakValue = signal[i];

        //peak is strong enough to possibly be a heartbeat
            if (peakValue > threshold) {
            //check if this peak could actually be a T-wave
                bool isTWave = false;

                if (lastPeak != -1 && i - lastPeak <= tWaveWindow) {
                    int slopeWindow = 0.075 * samplingRate;
                    double currentSlope = 0.0;
                    double previousSlope = 0.0;

                    if (i - slopeWindow >= 1 && lastPeak - slopeWindow >= 1) {
                        for (int j = i - slopeWindow + 1; j <= i; j++) {
                            currentSlope += abs(signal[j] - signal[j - 1]);
                    }

                    for (int j = lastPeak - slopeWindow + 1; j <= lastPeak; j++) {
                        previousSlope += abs(signal[j] - signal[j - 1]);
                    }

                    currentSlope /= slopeWindow;
                    previousSlope /= slopeWindow;

                    if (currentSlope < 0.5 * previousSlope) {
                        isTWave = true;
                    }
                }
            }

            //only accept it if it is not a T-wave
            if (!isTWave && (lastPeak == -1 || i - lastPeak > refractoryPeriod)) {
                signalPeakEstimate = 0.125 * peakValue + 0.875 * signalPeakEstimate;
                peaks.push_back(i);

                //calculate RR interval only if there was a previous heartbeat
                if (lastPeak != -1) {
                    int rr = i - lastPeak;
                    rrIntervals.push_back(rr);

                    if (rrIntervals.size() > 8) {
                        rrIntervals.erase(rrIntervals.begin());
                    }

                    double rrSum = 0.0;
                    for (int value : rrIntervals) {
                        rrSum += value;
                    }

                    rrAverage = rrSum / rrIntervals.size();
                    missedBeatLimit = 1.66 * rrAverage;
                }

                //current heartbeat becomes the new previous heartbeat
                lastPeak = i;
            }

        } else {
            //if peak too small = noise
            noisePeakEstimate = 0.125 * peakValue + 0.875 * noisePeakEstimate;

            //saving the strongest weak peak in case we missed a heartbeat
            if (peakValue > searchBackThreshold && peakValue > searchBackValue) {
                searchBackValue = peakValue;
                searchBackPeak = i;
            }
        }

        //update threshold for the next candidate peak
        threshold = noisePeakEstimate + 0.25 * (signalPeakEstimate - noisePeakEstimate);
        searchBackThreshold = 0.5 * threshold;
    }

    //search-back logic
    if (lastPeak != -1 && missedBeatLimit > 0 && i - lastPeak > missedBeatLimit) {
        if (searchBackPeak != -1) {
            int rr = searchBackPeak - lastPeak; //calculating RR interval for the recovered heartbeat
            rrIntervals.push_back(rr);

            if (rrIntervals.size() > 8) {
                rrIntervals.erase(rrIntervals.begin());
            }

            //recalculate the average RR interval
            double rrSum = 0.0;
            for (int value : rrIntervals) {
                rrSum += value;
            }

            rrAverage = rrSum / rrIntervals.size();
            missedBeatLimit = 1.66 * rrAverage;

            //add the recovered heartbeat
            peaks.push_back(searchBackPeak);
            lastPeak = searchBackPeak;

            searchBackPeak = -1;
            searchBackValue = 0.0;
             }
         }
    }

    return peaks;
}

