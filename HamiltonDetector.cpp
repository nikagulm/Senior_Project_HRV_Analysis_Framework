#include "HamiltonDetector.h"
#include <algorithm>
#include <cmath>

using namespace std;
HamiltonDetector::HamiltonDetector(int samplingRate) {
    this->samplingRate = samplingRate;
}


/*runs all Hamilton detector stages and returns the detected R-peaks
  @param ecgSignal original ECG samples.
  @return detected R-peak positions. */
vector<int> HamiltonDetector::detect(const vector<double>& ecgSignal) {
    if (ecgSignal.empty()) {
        return {};
    }

    //process the ECG through the HAMILTON STAGES - 4 main stages
    vector<double> filteredSignal = bandPassFilter(ecgSignal);
    vector<double> derivativeSignal = derivative(filteredSignal);
    vector<double> absoluteValues = absoluteSignal(derivativeSignal);
    vector<double> averagedSignal = movingAverage(absoluteValues);

    return detectPeaks(averagedSignal, ecgSignal);
}




/*  STAGE I : applies a simple band-pass filtering stage
   takes the original ECG signal and returns filtered ECG signal */
vector<double> HamiltonDetector::bandPassFilter(const vector<double>& signal) {
    vector<double> lowPass(signal.size(), 0.0);
    vector<double> highPass(signal.size(), 0.0);

    //calculate filter windows based on the sampling rate
    int lowPassWindow = static_cast<int>(0.05 * samplingRate);
    int highPassWindow = static_cast<int>(0.20 * samplingRate);

    if (lowPassWindow < 1) {
        lowPassWindow = 1;
    }

    if (highPassWindow < 1) {
        highPassWindow = 1;
    }

    double lowSum = 0.0;

    //apply the low-pass moving average
    for (int i = 0; i < signal.size(); i++) {
        lowSum += signal[i];
        if (i >= lowPassWindow) {
            lowSum -= signal[i - lowPassWindow];
        }

        lowPass[i] = lowSum / min(i + 1, lowPassWindow);
    }

    double highSum = 0.0;

    //remove the slower baseline component
    for (int i = 0; i < lowPass.size(); i++) {
        highSum += lowPass[i];
        if (i >= highPassWindow) {
            highSum -= lowPass[i - highPassWindow];
        }
        double baseline = highSum / min(i + 1, highPassWindow);
        highPass[i] = lowPass[i] - baseline;
    }
    return highPass;
}




/* STAGE 2: calculates the derivative of the filtered signal
   @param signal filtered ECG signal
   @return derivative values  */
vector<double> HamiltonDetector::derivative(const vector<double>& signal) {
    vector<double> result(signal.size(), 0.0);

    //measure the change between neighboring samples
    for (int i = 1; i < signal.size(); i++) {
        result[i] = signal[i] - signal[i - 1];
    }
    return result;
}




/* STAGE 3: converts all values to their absolute magnitude
   @param signal derivative signal
   @return absolute signal */
vector<double> HamiltonDetector::absoluteSignal(const vector<double>& signal) {
    vector<double> result(signal.size());
    for (int i = 0; i < signal.size(); i++) {
        result[i] = abs(signal[i]);
    }
    return result;
}



/* STAGE 4: smooths the signal using a moving average
   @param signal absolute derivative signal
   @return smoothed signal */
vector<double> HamiltonDetector::movingAverage(const vector<double>& signal) {
    //use an 80 ms moving-average window
    int windowSize = static_cast<int>(0.08 * samplingRate);

    if (windowSize < 1) {
        windowSize = 1;
    }
    vector<double> result(signal.size(), 0.0);
    double sum = 0.0;

    //use a running sum instead of recalculating the whole window
    for (int i = 0; i < signal.size(); i++) {
        sum += signal[i];

        if (i >= windowSize) {
            sum -= signal[i - windowSize];
        }

        if (i >= windowSize - 1) {
            result[i] = sum / windowSize;
        }
    }

    return result;
}




/* finds the actual R-peak close to a detected candidate
   @param signal original ECG signal
   @param candidatePosition detected candidate position
   @return R-peak position */
int HamiltonDetector::findRPeak(const vector<double>& signal, int candidatePosition) {
    //search around the candidate position for the actual R-peak
    int searchWindow = static_cast<int>(0.10 * samplingRate);

    int start = max(0, candidatePosition - searchWindow);
    int end = min(candidatePosition + searchWindow, static_cast<int>(signal.size()) - 1);

    if (start >= end) {
        return candidatePosition;
    }

    //calculate the local mean of the ECG
    double sum = 0.0;

    for (int i = start; i <= end; i++) {
        sum += signal[i];
    }

    double mean = sum / (end - start + 1);
    int rPeak = start;
    double largestDeviation = abs(signal[start] - mean);

    //choose the point furthest from the local mean
    for (int i = start + 1; i <= end; i++) {
        double deviation = abs(signal[i] - mean);

        if (deviation > largestDeviation) {
            largestDeviation = deviation;
            rPeak = i;
        }
    }

    return rPeak;
}




/* detects heartbeat candidates using an adaptive threshold
   @param processedSignal processed Hamilton signal
   @param originalSignal original ECG signal
   @return detected R-peak positions */
vector<int> HamiltonDetector::detectPeaks(const vector<double>& processedSignal, const vector<double>& originalSignal) {
    vector<int> peaks;
    if (processedSignal.empty()) {
        return peaks;
    }

    int initializationSamples = min(static_cast<int>(2 * samplingRate), static_cast<int>(processedSignal.size()));

    vector<double> signalPeaks;
    vector<double> noisePeaks;

    double maximumValue = 0.0;
    double sumValue = 0.0;

    //use the first two seconds to initialize the levels
    for (int i = 0; i < initializationSamples; i++) {
        sumValue += processedSignal[i];
        if (processedSignal[i] > maximumValue) {
            maximumValue = processedSignal[i];
        }
    }

    double averageValue = sumValue / initializationSamples;

    signalPeaks.push_back(maximumValue);
    noisePeaks.push_back(averageValue);
    double signalLevel = maximumValue;
    double noiseLevel = averageValue;

    double threshold = noiseLevel + 0.3 * (signalLevel - noiseLevel);

    int refractoryPeriod = static_cast<int>(0.20 * samplingRate);
    int tWaveWindow = static_cast<int>(0.36 * samplingRate);
    int slopeWindow = static_cast<int>(0.075 * samplingRate);

    int lastPeak = -1;

    for (int i = 1; i < processedSignal.size() - 1; i++) {
        //only check local maxima
        if (!(processedSignal[i] > processedSignal[i - 1] &&
              processedSignal[i] >= processedSignal[i + 1])) {
            continue;
        }

        double peakValue = processedSignal[i];
        if (peakValue > threshold) {
            bool isTWave = false;

            //compare slopes if the candidate appears soon after the last beat
            if (lastPeak != -1 &&
                i - lastPeak <= tWaveWindow &&
                i - slopeWindow >= 1 &&
                lastPeak - slopeWindow >= 1) {

                double currentSlope = 0.0;
                double previousSlope = 0.0;

                for (int j = i - slopeWindow + 1; j <= i; j++) {
                    currentSlope += abs(processedSignal[j] - processedSignal[j - 1]);
                }

                for (int j = lastPeak - slopeWindow + 1; j <= lastPeak; j++) {
                    previousSlope += abs(processedSignal[j] - processedSignal[j - 1]);
                }

                currentSlope /= slopeWindow;
                previousSlope /= slopeWindow;
                if (currentSlope < 0.5 * previousSlope) {
                    isTWave = true;
                }
            }

            if (!isTWave && (lastPeak == -1 || i - lastPeak > refractoryPeriod)) {
                int rPeak = findRPeak(originalSignal, i);

                if (peaks.empty() || rPeak - peaks.back() > refractoryPeriod) {
                    peaks.push_back(rPeak);
                }

                signalPeaks.push_back(peakValue);

                if (signalPeaks.size() > 8) {
                    signalPeaks.erase(signalPeaks.begin());
                }

                double signalSum = 0.0;
                for (double value : signalPeaks) {
                    signalSum += value;
                }

                signalLevel = signalSum / signalPeaks.size();
                lastPeak = i;
            }
        } else {
            noisePeaks.push_back(peakValue);

            if (noisePeaks.size() > 8) {
                noisePeaks.erase(noisePeaks.begin());
            }

            double noiseSum = 0.0;

            for (double value : noisePeaks) {
                noiseSum += value;
            }
            noiseLevel = noiseSum/noisePeaks.size();
        }

        //update the adaptive threshold
        threshold = noiseLevel + 0.3 * (signalLevel - noiseLevel);
    }

    return peaks;
}