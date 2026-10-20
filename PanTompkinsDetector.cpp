#include "PanTompkinsDetector.h"
#include <algorithm>
#include <cmath>

using namespace std;

/* creates a Pan-Tompkins detector.
   @param samplingRate original ECG sampling frequency in Hz.*/
PanTompkinsDetector::PanTompkinsDetector(int samplingRate) {
    originalSamplingRate = samplingRate;
}

/* runs all Pan-Tompkins stages and returns the detected R-peaks.
   @param ecgSignal original ECG samples.
   @return detected R-peak positions in original sample coordinates.*/
vector<int> PanTompkinsDetector::detect(const vector<double>& ecgSignal) {
    if (ecgSignal.empty()) {
        return {};
    }

    vector<double> processingSignal;

    //process the signal at 200 Hz
    if (originalSamplingRate != processingRate) {
        processingSignal = resample(ecgSignal, originalSamplingRate, processingRate);
    } else {
        processingSignal = ecgSignal;
    }

    vector<double> filteredSignal = bandPassFilter(processingSignal);
    vector<double> derivativeSignal = derivative(filteredSignal);
    vector<double> squaredSignal = squareSignal(derivativeSignal);
    vector<double> integratedSignal = movingWindowIntegration(squaredSignal);

    vector<int> processingPeaks = detectPeaks(integratedSignal, processingSignal);
    vector<int> originalPeaks;

    //convert peak positions back to the original sampling rate
    for (int peak : processingPeaks) {
        int convertedPeak = static_cast<int>(round(static_cast<double>(peak) * originalSamplingRate / processingRate));
        originalPeaks.push_back(convertedPeak);
    }

    return originalPeaks;
}



/* resamples the ECG signal using linear interpolation.
   @param signal ECG samples.
   @param originalRate original sampling frequency.
   @param targetRate target sampling frequency.
   @return resampled ECG signal. */
vector<double> PanTompkinsDetector::resample(const vector<double>& signal, int originalRate, int targetRate) {
    if (signal.empty() || originalRate <= 0 || targetRate <= 0) {
        return {};
    }

    int newSize = static_cast<int>(static_cast<double>(signal.size()) * targetRate / originalRate);
    vector<double> result(newSize);

    for (int i = 0; i < newSize; i++) {
        double originalPosition = static_cast<double>(i) * originalRate / targetRate;

        int left = static_cast<int>(originalPosition);
        int right = left + 1;

        if (right >= signal.size()) {
            right = signal.size() - 1;
        }

        double fraction = originalPosition - left;
        result[i] = signal[left] * (1.0 - fraction) + signal[right] * fraction;
    }

    return result;
}





/* applies the low-pass and high-pass filters used by Pan-Tompkins.
   @param signal ECG signal sampled at 200 Hz.
   @return filtered ECG signal. */
vector<double> PanTompkinsDetector::bandPassFilter(const vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }

    vector<double> lowPassValues(signal.size(), 0.0);

    //apply low-pass filter
    for (int i = 12; i < signal.size(); i++) {
        lowPassValues[i] = 2 * lowPassValues[i - 1] - lowPassValues[i - 2]
            + signal[i] - 2 * signal[i - 6] + signal[i - 12];
    }

    vector<double> highPassValues(signal.size(), 0.0);

    //apply high-pass filter
    for (int i = 32; i < signal.size(); i++) {
        highPassValues[i] = highPassValues[i - 1] - lowPassValues[i] / 32.0
            + lowPassValues[i - 16] - lowPassValues[i - 17]
            + lowPassValues[i - 32] / 32.0;
    }

    return highPassValues;
}



/*calculates the derivative to emphasize rapid changes in the ECG.
  @param signal filtered ECG signal.
  @return derivative values. */
vector<double> PanTompkinsDetector::derivative(const vector<double>& signal) {
    if (signal.size() < 5) {
        return {};
    }

    vector<double> derivativeValues(signal.size(), 0.0);

    for (int i = 2; i < signal.size() - 2; i++) {
        derivativeValues[i] = (-signal[i - 2] - 2 * signal[i - 1]
            + 2 * signal[i + 1] + signal[i + 2]) / 8.0;
    }

    return derivativeValues;
}




/*squares the derivative values to make large changes more noticeable.
  @param signal derivative signal.
  @return squared signal. */
vector<double> PanTompkinsDetector::squareSignal(const vector<double>& signal) {
    vector<double> squaredValues(signal.size());

    for (int i = 0; i < signal.size(); i++) {
        squaredValues[i] = signal[i] * signal[i];
    }

    return squaredValues;
}



/* applies a 150 ms moving-window integration.
   @param signal squared ECG signal.
   @return integrated signal. */
vector<double> PanTompkinsDetector::movingWindowIntegration(const vector<double>& signal) {
    if (signal.empty()) {
        return {};
    }

    int windowSize = static_cast<int>(0.15 * processingRate);
    vector<double> integratedValues(signal.size(), 0.0);
    double sum = 0.0;
    for (int i = 0; i < signal.size(); i++) {
        sum += signal[i];

        if (i >= windowSize) {
            sum -= signal[i - windowSize];
        }

        if (i >= windowSize - 1) {
            integratedValues[i] = sum / windowSize;
        }
    }
    return integratedValues;
}





/* finds the actual R-peak close to a detected QRS candidate.
   - the local mean is removed so the ECG baseline does not affect which point is selected as the R-peak.
   @param signal raw ECG signal sampled at 200 Hz.
   @param candidatePosition position found in the integrated signal.
   @return position of the R-peak.*/
int PanTompkinsDetector::findRPeak(const vector<double>& signal, int candidatePosition) {
    int searchWindow = static_cast<int>(0.20 * processingRate);

    int start = max(0, candidatePosition - searchWindow);
    int end = min(candidatePosition, static_cast<int>(signal.size()) - 1);

    if (start >= end) {
        return candidatePosition;
    }

    //calculate the local baseline
    double sum = 0.0;

    for (int i = start; i <= end; i++) {
        sum += signal[i];
    }

    double mean = sum / (end - start + 1);

    int rPeak = start;
    double largestDeviation = abs(signal[start] - mean);

    //find the point furthest from the local baseline
    for (int i = start + 1; i <= end; i++) {
        double deviation = abs(signal[i] - mean);

        if (deviation > largestDeviation) {
            largestDeviation = deviation;
            rPeak = i;
        }
    }

    return rPeak;
}




/* detects QRS candidates using adaptive thresholds.
   it also handles the refractory period, possible T-waves and
   search-back for beats that may have been missed.
 
 -> @param integratedSignal moving-window integrated signal.
 -> @param originalSignal raw ECG signal sampled at 200 Hz.
 -> @return detected R-peak positions.  */
vector<int> PanTompkinsDetector::detectPeaks(const vector<double>& integratedSignal, const vector<double>& originalSignal) {
    vector<int> peaks;

    if (integratedSignal.empty() || originalSignal.empty()) {
        return peaks;
    }

    //use the first two seconds to initialize the thresholds
    int initializationSamples = 2 * processingRate;

    if (initializationSamples > integratedSignal.size()) {
        initializationSamples = integratedSignal.size();
    }

    double maximumValue = 0.0;
    double sumValue = 0.0;

    for (int i = 0; i < initializationSamples; i++) {
        sumValue += integratedSignal[i];

        if (integratedSignal[i] > maximumValue) {
            maximumValue = integratedSignal[i];
        }
    }

    double averageValue = sumValue / initializationSamples;
    double signalPeakEstimate = 0.25 * maximumValue;
    double noisePeakEstimate = 0.5 * averageValue;

    double threshold = noisePeakEstimate + 0.25 * (signalPeakEstimate - noisePeakEstimate);
    double searchBackThreshold = 0.5 * threshold;

    int refractoryPeriod = static_cast<int>(0.20 * processingRate);
    int tWaveWindow = static_cast<int>(0.36 * processingRate);
    int slopeWindow = static_cast<int>(0.075 * processingRate);

    int lastCandidate = -1;
    int lastRPeak = -1;

    vector<int> rrIntervals;

    double rrAverage = 0.0;
    double missedBeatLimit = 0.0;

    int searchBackPeak = -1;
    double searchBackValue = 0.0;

    //check local maxima and decide whether they are QRS complexes
    for (int i = 1; i < integratedSignal.size() - 1; i++) {
        if (!(integratedSignal[i] > integratedSignal[i - 1] &&
              integratedSignal[i] >= integratedSignal[i + 1])) {
            continue;
        }

        double peakValue = integratedSignal[i];

        if (peakValue > threshold) {
            bool isTWave = false;

            //check if the candidate could be a T-wave
            if (lastCandidate != -1 &&
                i - lastCandidate <= tWaveWindow &&
                i - slopeWindow >= 1 &&
                lastCandidate - slopeWindow >= 1) {

                double currentSlope = 0.0;
                double previousSlope = 0.0;

                for (int j = i - slopeWindow + 1; j <= i; j++) {
                    currentSlope += abs(integratedSignal[j] - integratedSignal[j - 1]);
                }

                for (int j = lastCandidate - slopeWindow + 1; j <= lastCandidate; j++) {
                    previousSlope += abs(integratedSignal[j] - integratedSignal[j - 1]);
                }

                currentSlope /= slopeWindow;
                previousSlope /= slopeWindow;

                if (currentSlope < 0.5 * previousSlope) {
                    isTWave = true;
                }
            }

            if (!isTWave && (lastCandidate == -1 || i - lastCandidate > refractoryPeriod)) {
                int rPeak = findRPeak(originalSignal, i);

                if (lastRPeak == -1 || rPeak - lastRPeak > refractoryPeriod) {
                    peaks.push_back(rPeak);
                    lastRPeak = rPeak;
                }

                signalPeakEstimate = 0.125 * peakValue + 0.875 * signalPeakEstimate;

                //update RR interval information
                if (lastCandidate != -1) {
                    int rr = i - lastCandidate;
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

                lastCandidate = i;
                searchBackPeak = -1;
                searchBackValue = 0.0;
            }
        } else {
            noisePeakEstimate = 0.125 * peakValue + 0.875 * noisePeakEstimate;

            if (peakValue > searchBackThreshold && peakValue > searchBackValue) {
                searchBackValue = peakValue;
                searchBackPeak = i;
            }
        }

        //update the adaptive thresholds
        threshold = noisePeakEstimate + 0.25 * (signalPeakEstimate - noisePeakEstimate);
        searchBackThreshold = 0.5 * threshold;

        //search back when a heartbeat may have been missed
        if (lastCandidate != -1 &&
            missedBeatLimit > 0 &&
            i - lastCandidate > missedBeatLimit &&
            searchBackPeak != -1 &&
            searchBackPeak > lastCandidate + refractoryPeriod) {

            int rPeak = findRPeak(originalSignal, searchBackPeak);

            if (lastRPeak == -1 || rPeak - lastRPeak > refractoryPeriod) {
                peaks.push_back(rPeak);
                lastRPeak = rPeak;
            }

            int rr = searchBackPeak - lastCandidate;
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

            lastCandidate = searchBackPeak;
            signalPeakEstimate = 0.25 * searchBackValue + 0.75 * signalPeakEstimate;

            searchBackPeak = -1;
            searchBackValue = 0.0;

            threshold = noisePeakEstimate + 0.25 * (signalPeakEstimate - noisePeakEstimate);
            searchBackThreshold = 0.5 * threshold;
        }
    }

    sort(peaks.begin(), peaks.end());

    return peaks;
}