#include "HamiltonDetector.h"
#include <algorithm>
#include <cmath>
using namespace std;

HamiltonDetector::HamiltonDetector(int samplingRate, bool refineRPeak) {
    originalSamplingRate = samplingRate;
    this->refineRPeak = refineRPeak;
}


//runs all Hamilton stages and returns the detected R-peaks
vector<int> HamiltonDetector::detect(const vector<double>& ecgSignal) {
    if (ecgSignal.empty()) {
        return {};
    }

    //STEP 1: change to 200 Hz, the rate the filters are designed for
    vector<double> signal = ecgSignal;
    if (originalSamplingRate != processingRate) {
        signal = resample(ecgSignal, originalSamplingRate, processingRate);
    }


    //Hamilton's code expects ADC units (200 units = 1 mV)
    vector<double> adcSignal(signal.size());
    for (int i = 0; i < signal.size(); i++) {
        adcSignal[i] = signal[i] * 200.0;
    }

    //STEP 2: band-pass filter = low-pass, then high-pass
    vector<double> lowPassed = lowPassFilter(adcSignal);
    vector<double> bandPassed = highPassFilter(lowPassed);

    //STEP 3: derivative, absolute value and moving-window integration
    vector<double> derivativeSignal = derivative(bandPassed);
    vector<double> absoluteValues = absoluteSignal(derivativeSignal);
    vector<double> integratedSignal = movingWindowIntegration(absoluteValues);

    //slope of the unfiltered ECG (needed for the baseline-shift and T-wave rules)
    vector<double> rawSlope = derivative(adcSignal);

    //STEP 4: decide which peaks are heartbeats
    vector<int> qrs = findQRSComplexes(integratedSignal, rawSlope);

    //STEP 5: convert the positions back to the original sampling rate
    vector<int> rPeaks;

    for (int position : qrs) {
        int rPeak = max(0, position);

        //optional: move to the exact R-peak
        if (refineRPeak) {
            rPeak = findRPeak(signal, position);
        }

        int originalPosition = static_cast<int>(round(static_cast<double>(rPeak) * originalSamplingRate / processingRate));

        //avoid adding the same R-peak twice
        if (rPeaks.empty() || originalPosition > rPeaks.back()) {
            rPeaks.push_back(originalPosition);
        }
    }

    return rPeaks;
}



/* changes the sampling rate (360 Hz -> 200 Hz).
   each new sample is estimated from the two original samples around it (linear interpolation) */
vector<double> HamiltonDetector::resample(const vector<double>& signal, int originalRate, int targetRate) {
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



/* returns a sample, treating samples before the start of the signal as 0 */
double HamiltonDetector::sampleAt(const vector<double>& v, int i) {
    if (i < 0) {
        return 0.0;
    }
    return v[i];
}



/* low-pass filter, the same as in Pan-Tompkins (delay 5 samples):
   y(n) = 2y(n-1) - y(n-2) + x(n) - 2x(n-6) + x(n-12)
   the output is divided by the filter gain (36) */
vector<double> HamiltonDetector::lowPassFilter(const vector<double>& x) {
    vector<double> y(x.size(), 0.0);
    vector<double> output(x.size(), 0.0);

    for (int n = 0; n < x.size(); n++) {
        y[n] = 2 * sampleAt(y, n - 1) - sampleAt(y, n - 2)
             + x[n] - 2 * sampleAt(x, n - 6) + sampleAt(x, n - 12);
        output[n] = y[n] / 36.0;
    }

    return output;
}



/* high-pass filter (delay 13 samples): subtracts the 26-sample (125 ms)
   moving average from the signal, which removes baseline wander
       sum(n) = sum(n-1) + x(n) - x(n-26)
       z(n)   = x(n-13) - sum(n) / 26 */
vector<double> HamiltonDetector::highPassFilter(const vector<double>& x) {
    int length = 26;
    vector<double> z(x.size(), 0.0);
    double sum = 0.0;

    for (int n = 0; n < x.size(); n++) {
        sum += x[n] - sampleAt(x, n - length);
        z[n] = sampleAt(x, n - length / 2) - sum / length;
    }

    return z;
}



/* derivative over 3 samples (about 15 ms): y(n) = x(n) - x(n-3) */
vector<double> HamiltonDetector::derivative(const vector<double>& x) {
    vector<double> y(x.size(), 0.0);

    for (int n = 0; n < x.size(); n++) {
        y[n] = x[n] - sampleAt(x, n - 3);
    }

    return y;
}



/* takes the absolute value of every sample.
   (Hamilton uses this instead of the squaring used in Pan-Tompkins) */
vector<double> HamiltonDetector::absoluteSignal(const vector<double>& x) {
    vector<double> y(x.size());

    for (int n = 0; n < x.size(); n++) {
        y[n] = abs(x[n]);
    }

    return y;
}



/* moving-window integration over 80 ms (17 samples)
   as in Hamilton's code, the output is limited to 32000 */
vector<double> HamiltonDetector::movingWindowIntegration(const vector<double>& x) {
    vector<double> y(x.size(), 0.0);
    double sum = 0.0;

    for (int n = 0; n < x.size(); n++) {
        //add the new sample and remove the one that left the window
        sum += x[n] - sampleAt(x, n - windowWidth);
        y[n] = min(sum / windowWidth, 32000.0);
    }

    return y;
}




//DECISION RULES
/* goes through the integrated signal one sample at a time and decides which peaks are 
heartbeats, using the five rules:
     1. ignore peaks within 200 ms of a larger peak
     2. ignore baseline shifts
     3. a weak peak within 360 ms of a beat is a T-wave
     4. a peak above the detection threshold is a QRS, otherwise noise
     5. search back for a missed beat after 1.5 x the RR interval */
vector<int> HamiltonDetector::findQRSComplexes(const vector<double>& integrated, const vector<double>& rawSlope) {
    //reset everything in case detect() is called more than once
    peakMax = 0.0;
    timeSinceMax = 0;
    lastValue = 0.0;
    heldPeak = 0.0;
    blankCounter = 0;
    qrsPeaks.assign(8, 0.0);
    noisePeaks.assign(8, 0.0);
    rrIntervals.assign(8, ms1000);
    qrsMedian = 0.0;
    noiseMedian = 0.0;
    rrMedian = ms1000;
    detectionThreshold = 0.0;
    count = 0;
    searchBackLimit = ms1500;
    searchBackPeak = 0.0;
    searchBackLocation = 0;
    searchBackSlope = 0.0;
    lastQRSSlope = 0.0;
    learnedPeaks = 0;
    secondCounter = 0;
    secondMax = 0.0;
    resetPeaks.clear();
    qrsPositions.clear();

    for (int n = 0; n < integrated.size(); n++) {

        //find peaks
        double peak = findPeak(integrated[n]);

        //rule 1 keep only the largest peak within 200 ms
        double newPeak = 0.0;

        if (peak > 0 && blankCounter == 0) {
            //start holding a new peak
            heldPeak = peak;
            blankCounter = ms200;
        } else if (peak > 0 && peak > heldPeak) {
            //a bigger peak replaces the held one, and the 200 ms start again
            heldPeak = peak;
            blankCounter = ms200;
        } else if (blankCounter > 0) {
            blankCounter--;

            //no bigger peak came within 200 ms, so the held peak is final
            if (blankCounter == 0) {
                newPeak = heldPeak;
            }
        }

        count++;

        if (learnedPeaks < 8) {
            //learning phase first 8 seconds
            if (newPeak > 0) {
                count = windowWidth;
            }

            secondCounter++;

            //every second save the largest peak of that second as a QRS peak
            if (secondCounter == ms1000) {
                qrsPeaks[learnedPeaks] = secondMax;
                learnedPeaks++;
                secondCounter = 0;
                secondMax = 0.0;

                //after 8 seconds set the first threshold
                if (learnedPeaks == 8) {
                    qrsMedian = median(qrsPeaks);
                    noiseMedian = 0.0;
                    rrMedian = ms1000;
                    searchBackLimit = ms1500 + ms150;
                    updateThreshold();
                }
            }

            if (newPeak > secondMax) {
                secondMax = newPeak;
            }
        } else {
            //detection phase
            double maxSlope = 0.0;
            bool baselineShift = false;
            bool tWave = false;

            if (newPeak > 0) {
                //rule 2 a real QRS has both a positive and a negative slope
                baselineShift = isBaselineShift(rawSlope, n, maxSlope);

                //rule 3 within 360 ms of the last QRS, a peak with less than
                //half of that QRS's slope is a T-wave
                tWave = (count - windowWidth < ms360) && (maxSlope < 0.5 * lastQRSSlope);
            }

            //a baseline shift is ignored completely
            if (newPeak > 0 && !baselineShift) {

                //rule 4 above the threshold (and not a T-wave) = QRS
                if (newPeak > detectionThreshold && !tWave) {
                    acceptQRS(newPeak, count - windowWidth);

                    //the QRS happened "windowWidth + filterDelay" samples ago
                    qrsPositions.push_back(n - (windowWidth + filterDelay));
                    lastQRSSlope = maxSlope;
                    count = windowWidth;
                } else {
                    //otherwise noise (a T-wave also counts as noise)
                    noisePeaks.insert(noisePeaks.begin(), newPeak);
                    noisePeaks.pop_back();
                    noiseMedian = median(noisePeaks);
                    updateThreshold();

                    //remember the largest noise peak for search-back,
                    //but only if it is at least 360 ms after the last QRS
                    if (newPeak > searchBackPeak && count - windowWidth >= ms360) {
                        searchBackPeak = newPeak;
                        searchBackLocation = count - windowWidth;
                        searchBackSlope = maxSlope;
                    }
                }
            }

            //rule 5 search-back 
            //no QRS for 1.5 x RR interval. accept the saved peak if it is above half the threshold
            if (useSearchBack() && count > searchBackLimit && searchBackPeak > 0.5 * detectionThreshold) {
                acceptQRS(searchBackPeak, searchBackLocation);
                count = count - searchBackLocation;
                qrsPositions.push_back(n - (count + filterDelay));
                lastQRSSlope = searchBackSlope;
            }
        }

        //reset: no QRS for 8 seconds -> learn the levels again 
        if (learnedPeaks == 8) {
            secondCounter++;

            if (secondCounter == ms1000) {
                resetPeaks.push_back(secondMax);
                secondCounter = 0;
                secondMax = 0.0;

                if (resetPeaks.size() == 8) {
                    qrsPeaks = resetPeaks;
                    noisePeaks.assign(8, 0.0);
                    qrsMedian = median(resetPeaks);
                    noiseMedian = 0.0;
                    rrMedian = ms1000;
                    searchBackLimit = ms1500 + ms150;
                    updateThreshold();
                    restartLearning();
                    searchBackPeak = 0.0;
                }
            }

            if (newPeak > secondMax) {
                secondMax = newPeak;
            }
        }
    }

    sort(qrsPositions.begin(), qrsPositions.end());
    return qrsPositions;
}



/* finds peaks in the integrated signal (called once per sample)
   a peak is reported when the signal drops below half of its maximum
   or when 95 ms have passed since the maximum. */
double HamiltonDetector::findPeak(double value) {
    double peak = 0.0;

    if (timeSinceMax > 0) {
        timeSinceMax++;
    }

    if (value > lastValue && value > peakMax) {
        //the signal is still rising
        peakMax = value;
        if (peakMax > 2) {
            timeSinceMax = 1;
        }
    } else if (value < peakMax / 2) {
        //the signal dropped below half of the maximum
        peak = peakMax;
        peakMax = 0.0;
        timeSinceMax = 0;
    } else if (timeSinceMax > ms95) {
        //95 ms passed since the maximum
        peak = peakMax;
        peakMax = 0.0;
        timeSinceMax = 0;
    }

    lastValue = value;
    return peak;
}







/* rule 2 checks if a peak is a baseline shift instead of a QRS
   looks at the slope of the unfiltered ECG in a 220 ms window around the peak
   a real QRS has a large positive AND a large negative slope:
   - each must be at least 1/8 of the other
   - they must be less than 150 ms apart
   also saves the largest positive slope in maxSlope (used by the T-wave rule). */
bool HamiltonDetector::isBaselineShift(const vector<double>& rawSlope, int n, double& maxSlope) {
    //the QRS lies in the oldest 220 ms of the last 96 samples
    int start = n - (slopeBufferLength - 1);

    double maxValue = 0.0;
    double minValue = 0.0;
    int maxTime = 0;
    int minTime = 0;

    //find the largest positive and the largest negative slope
    for (int t = 0; t < ms220; t++) {
        double value = sampleAt(rawSlope, start + t);

        if (value > maxValue) {
            maxValue = value;
            maxTime = t;
        } else if (value < minValue) {
            minValue = value;
            minTime = t;
        }
    }

    maxSlope = maxValue;
    minValue = -minValue;

    bool bothSlopes = maxValue > minValue / 8 && minValue > maxValue / 8;
    bool closeTogether = abs(maxTime - minTime) < ms150;

    return !(bothSlopes && closeTogether);
}



/* median of the values, as in Hamilton's code (for 8 values: the 5th smallest). */
double HamiltonDetector::median(const vector<double>& values) {
    vector<double> sorted = values;
    sort(sorted.begin(), sorted.end());
    return sorted[sorted.size() / 2];
}



/* recalculates the detection threshold
   DT = noise level + 0.475 * (QRS level - noise level) */
void HamiltonDetector::updateThreshold() {
    detectionThreshold = noiseMedian + thresholdCoefficient * (qrsMedian - noiseMedian);
}



/* updates the QRS level, the RR median and the search-back limit after a QRS */
void HamiltonDetector::acceptQRS(double peakHeight, int rrInterval) {
    //QRS level = median of the last 8 QRS peaks
    qrsPeaks.insert(qrsPeaks.begin(), peakHeight);
    qrsPeaks.pop_back();
    qrsMedian = median(qrsPeaks);
    updateThreshold();

    //RR median = median of the last 8 RR intervals
    rrIntervals.insert(rrIntervals.begin(), rrInterval);
    rrIntervals.pop_back();
    rrMedian = median(rrIntervals);

    //search back if no QRS comes within 1.5 x RR median
    searchBackLimit = static_cast<int>(rrMedian * 1.5) + windowWidth;

    searchBackPeak = 0.0;
    restartLearning();
}



// restarts the 8-second reset timer (called after every QRS) 
void HamiltonDetector::restartLearning() {
    secondCounter = 0;
    secondMax = 0.0;
    resetPeaks.clear();
}




/* finds the exact R-peak near a detection (my own addition)
     picks the point furthest from the local baseline */
int HamiltonDetector::findRPeak(const vector<double>& signal, int approximatePosition) {
    int size = signal.size();
    int searchWindow = static_cast<int>(0.075 * processingRate);    //75 ms each side
    int baselineWindow = static_cast<int>(0.20 * processingRate);   //200 ms each side

    int start = max(0, approximatePosition - searchWindow);
    int end = min(size - 1, approximatePosition + searchWindow);

    if (start >= end) {
        return max(0, min(size - 1, approximatePosition));
    }

    //calculate the local baseline
    int baselineStart = max(0, approximatePosition - baselineWindow);
    int baselineEnd = min(size - 1, approximatePosition + baselineWindow);

    double sum = 0.0;
    for (int i = baselineStart; i <= baselineEnd; i++) {
        sum += signal[i];
    }
    double mean = sum / (baselineEnd - baselineStart + 1);

    //find the point furthest from the local baseline
    int rPeak = start;
    for (int i = start + 1; i <= end; i++) {
        if (abs(signal[i] - mean) > abs(signal[rPeak] - mean)) {
            rPeak = i;
        }
    }

    return rPeak;
}