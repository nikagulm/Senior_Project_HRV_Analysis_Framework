#include "PanTompkinsDetector.h"
#include <algorithm>
#include <cmath>
using namespace std;

PanTompkinsDetector::PanTompkinsDetector(int samplingRate) {
    originalSamplingRate = samplingRate;
}


/* runs all Pan-Tompkins stages and returns the detected R-peaks
   takes ecgSignal original ECG samples
   returns detected R-peak positions in original sample coordinates */
vector<int> PanTompkinsDetector::detect(const vector<double>& ecgSignal) {
    if (ecgSignal.empty()) {
        return {};
    }

    //STEP 1: the published filters are designed for 200 Hz
    vector<double> signal = ecgSignal;
    if (originalSamplingRate != processingRate) {
        signal = resample(ecgSignal, originalSamplingRate, processingRate);
    }

    //STEP 2: band-pass filter = low-pass followed by high-pass
    vector<double> lowPassed = lowPassFilter(signal);
    vector<double> bandPassed = highPassFilter(lowPassed);

    //STEP 3: derivative, squaring and moving-window integration
    vector<double> derivativeSignal = derivative(bandPassed);
    vector<double> squaredSignal = squareSignal(derivativeSignal);
    vector<double> integratedSignal = movingWindowIntegration(squaredSignal);

    //STEP 4: decision rules
    vector<int> qrs = findQRSComplexes(integratedSignal, bandPassed);

    //STEP 5: locate each R-peak and convert back to the original sampling rate
    vector<int> rPeaks;

    for (int position : qrs) {
        //the band-passed signal is delayed by the two filters
        int approximatePosition = position - (lowPassDelay + highPassDelay);
        int rPeak = findRPeak(signal, approximatePosition);

        int originalPosition = static_cast<int>(round(static_cast<double>(rPeak) * originalSamplingRate / processingRate));

        //avoid adding the same R-peak twice
        if (rPeaks.empty() || originalPosition > rPeaks.back()) {
            rPeaks.push_back(originalPosition);
        }
    }

    return rPeaks;
}



/* changes the sampling rate from 360 Hz (dataset) to 200 Hz (Pan-Tompkins filters are designed for 200 Hz)
   each new sample is estimated from the two original samples around it (linear interpolation). */
vector<double> PanTompkinsDetector::resample(const vector<double>& signal, int originalRate, int targetRate) {
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



/* returns v[i], or 0 if i is before the start of the signal
   the filters below treat samples before the start as zero. */
double PanTompkinsDetector::sampleAt(const vector<double>& v, int i) {
    if (i < 0) {
        return 0.0;
    }
    return v[i];
}



/* low-pass filter from the paper (cutoff about 11 Hz, gain 36, delay 5 samples):
   y(n) = 2y(n-1) - y(n-2) + x(n) - 2x(n-6) + x(n-12) */
vector<double> PanTompkinsDetector::lowPassFilter(const vector<double>& x) {
    vector<double> y(x.size(), 0.0);

    for (int n = 0; n < x.size(); n++) {
        y[n] = 2 * sampleAt(y, n - 1) - sampleAt(y, n - 2)
             + x[n] - 2 * sampleAt(x, n - 6) + sampleAt(x, n - 12);
    }

    return y;
}



/* high-pass filter from the paper (cutoff about 5 Hz, gain 32, delay 16 samples).
   it is an all-pass (delay of 16 samples) minus a 32-sample low-pass:
       sum(n) = sum(n-1) + x(n) - x(n-32)
       y(n)   = 32x(n-16) - sum(n)  */
vector<double> PanTompkinsDetector::highPassFilter(const vector<double>& x) {
    vector<double> sum(x.size(), 0.0);
    vector<double> y(x.size(), 0.0);

    for (int n = 0; n < x.size(); n++) {
        sum[n] = sampleAt(sum, n - 1) + x[n] - sampleAt(x, n - 32);
        y[n] = 32 * sampleAt(x, n - 16) - sum[n];
    }

    return y;
}



/* five-point derivative from the paper (delay 2 samples):
   y(n) = (1/8T) * [2x(n) + x(n-1) - x(n-3) - 2x(n-4)],  T = 1/200 s */
vector<double> PanTompkinsDetector::derivative(const vector<double>& x) {
    vector<double> y(x.size(), 0.0);
    double T = 1.0 / processingRate;

    for (int n = 0; n < x.size(); n++) {
        y[n] = (2 * x[n] + sampleAt(x, n - 1) - sampleAt(x, n - 3) - 2 * sampleAt(x, n - 4)) / (8 * T);
    }

    return y;
}



/* squares every value so all values are positive and large slopes stand out */
vector<double> PanTompkinsDetector::squareSignal(const vector<double>& x) {
    vector<double> y(x.size());

    for (int n = 0; n < x.size(); n++) {
        y[n] = x[n] * x[n];
    }

    return y;
}



/* moving-window integration with a 150 ms window (30 samples at 200 Hz):
   y(n) = (1/N) * [x(n-(N-1)) + ... + x(n)] */
vector<double> PanTompkinsDetector::movingWindowIntegration(const vector<double>& x) {
    int windowSize = static_cast<int>(0.15 * processingRate);
    vector<double> y(x.size(), 0.0);

    double sum = 0.0;
    for (int n = 0; n < x.size(); n++) {
        sum += x[n];

        //remove the sample that left the window
        if (n >= windowSize) {
            sum -= x[n - windowSize];
        }

        y[n] = sum / windowSize;
    }

    return y; //integrated signal
}




//   -- DECISION RULES 


/* goes through all peaks of the integrated signal and decides
   whether each one is a QRS complex or noise.
   returns QRS positions in the band-passed signal. */
vector<int> PanTompkinsDetector::findQRSComplexes(const vector<double>& integrated, const vector<double>& bandPassed) {
    //reset everything in case detect() is called more than once
    qrsPositions.clear();
    noisePeaks.clear();
    recentRR.clear();
    regularRR.clear();
    rrAverage1 = 0.0;
    rrAverage2 = 0.0;
    rrMissedLimit = 0.0;
    irregularRhythm = false;
    lastQRSPosition = -1;
    lastQRSSlope = 0.0;

    //learning phase 1: the first 2 seconds only set the starting thresholds,
    //no QRS complexes are reported during this time
    int learningEnd = learnThresholds(integrated, bandPassed);

    int refractoryPeriod = static_cast<int>(0.20 * processingRate);   //200 ms

    for (int i = max(1, learningEnd); i < integrated.size() - 1; i++) {
        //a peak is a point higher than its neighbours
        bool isPeak = integrated[i] > integrated[i - 1] && integrated[i] >= integrated[i + 1];
        if (!isPeak) {
            continue;
        }

        //if no QRS was found for too long, look back for a missed beat
        if (lastQRSPosition != -1 && rrMissedLimit > 0) {
            while (i - lastQRSPosition > rrMissedLimit) {
                bool found = searchBack();
                if (!found) {
                    break;
                }
            }
        }

        //no new QRS can start within 200 ms of the previous one
        if (lastQRSPosition != -1 && i - lastQRSPosition <= refractoryPeriod) {
            continue;
        }

        PeakCandidate candidate = makeCandidate(i, integrated, bandPassed);

        //a QRS must pass the threshold in BOTH signals
        bool aboveThresholdI = candidate.peakI > thresholdI1;
        bool aboveThresholdF = candidate.peakF > thresholdF1;

        if (aboveThresholdI && aboveThresholdF) {
            if (isTWave(candidate)) {
                //a T-wave is treated like any other noise peak
                markAsNoise(candidate);
            } else {
                //signal peak update from the paper
                signalPeakI = 0.125 * candidate.peakI + 0.875 * signalPeakI;
                signalPeakF = 0.125 * candidate.peakF + 0.875 * signalPeakF;
                acceptQRS(candidate);
            }
        } else {
            markAsNoise(candidate);
        }

        updateThresholds();
    }

    sort(qrsPositions.begin(), qrsPositions.end());
    return qrsPositions;
}



/* collects the values the decision rules need for one peak
   the matching band-passed peak is searched in the window that formed
   the integrated peak (150 ms integration + 2 samples derivative delay)
   returns filled candidate. */
PeakCandidate PanTompkinsDetector::makeCandidate(int position, const vector<double>& integrated, const vector<double>& bandPassed) {
    PeakCandidate candidate;
    candidate.position = position;
    candidate.peakI = integrated[position];
    candidate.peakF = 0.0;
    candidate.positionF = position;
    candidate.slope = 0.0;

    int window = static_cast<int>(0.15 * processingRate) + 2;
    int start = max(1, position - window);

    for (int j = start; j <= position; j++) {
        //largest band-passed value (absolute, so negative QRS also work)
        if (abs(bandPassed[j]) > candidate.peakF) {
            candidate.peakF = abs(bandPassed[j]);
            candidate.positionF = j;
        }

        //steepest slope in the window
        double slope = abs(bandPassed[j] - bandPassed[j - 1]);
        if (slope > candidate.slope) {
            candidate.slope = slope;
        }
    }

    return candidate;
}



/* learning phase 1: sets the starting signal and noise levels
   from the peaks found in the first 2 seconds.
   - the highest peak is taken as the first signal peak (SPK)
   - the average of the other peaks is taken as the first noise peak (NPK)
   (the paper says the thresholds come from the signal and noise peaks found
   during learning, but does not give the exact formula, so this is my choice)
   return position where normal detection starts */
int PanTompkinsDetector::learnThresholds(const vector<double>& integrated, const vector<double>& bandPassed) {
    int end = min(static_cast<int>(integrated.size()) - 1, 2 * processingRate);

    vector<PeakCandidate> learningPeaks;
    for (int i = 1; i < end; i++) {
        bool isPeak = integrated[i] > integrated[i - 1] && integrated[i] >= integrated[i + 1];
        if (isPeak) {
            learningPeaks.push_back(makeCandidate(i, integrated, bandPassed));
        }
    }

    if (learningPeaks.empty()) {
        updateThresholds();
        return end;
    }

    //find the highest peak
    int highest = 0;
    for (int k = 1; k < learningPeaks.size(); k++) {
        if (learningPeaks[k].peakI > learningPeaks[highest].peakI) {
            highest = k;
        }
    }

    signalPeakI = learningPeaks[highest].peakI;
    signalPeakF = learningPeaks[highest].peakF;

    //average of all the other peaks
    double sumI = 0.0;
    double sumF = 0.0;
    int count = 0;
    for (int k = 0; k < learningPeaks.size(); k++) {
        if (k != highest) {
            sumI += learningPeaks[k].peakI;
            sumF += learningPeaks[k].peakF;
            count++;
        }
    }

    if (count > 0) {
        noisePeakI = sumI / count;
        noisePeakF = sumF / count;
    }

    updateThresholds();
    return end;
}



/* recalculates all four thresholds from the current signal and noise levels
   THRESHOLD1 = NPK + 0.25 * (SPK - NPK)
   THRESHOLD2 = 0.5 * THRESHOLD1
   if the rhythm is irregular, THRESHOLD1 is halved */
void PanTompkinsDetector::updateThresholds() {
    thresholdI1 = noisePeakI + 0.25 * (signalPeakI - noisePeakI);
    thresholdI2 = 0.5 * thresholdI1;

    thresholdF1 = noisePeakF + 0.25 * (signalPeakF - noisePeakF);
    thresholdF2 = 0.5 * thresholdF1;

    if (irregularRhythm) {
        thresholdI1 = 0.5 * thresholdI1;
        thresholdF1 = 0.5 * thresholdF1;
    }
}



/* updates the two RR averages after a new beat.
   RR AVERAGE1 = mean of the last 8 RR intervals
   RR AVERAGE2 = mean of the last 8 RR intervals between 92% and 116% of RR AVERAGE2
   RR MISSED LIMIT = 166% of RR AVERAGE2, or of RR AVERAGE1 after a sudden rate change */
void PanTompkinsDetector::updateRRAverages(int rrInterval) {
    //RR AVERAGE1
    recentRR.push_back(rrInterval);
    if (recentRR.size() > 8) {
        recentRR.erase(recentRR.begin());
    }

    double sum = 0.0;
    for (int value : recentRR) {
        sum += value;
    }
    rrAverage1 = sum / recentRR.size();

    //RR AVERAGE2
    if (regularRR.empty()) {
        //learning phase 2: the first RR interval starts the average
        regularRR.push_back(rrInterval);
    } else if (rrInterval >= 0.92 * rrAverage2 && rrInterval <= 1.16 * rrAverage2) {
        regularRR.push_back(rrInterval);
        if (regularRR.size() > 8) {
            regularRR.erase(regularRR.begin());
        }
    }

    sum = 0.0;
    for (int value : regularRR) {
        sum += value;
    }
    rrAverage2 = sum / regularRR.size();

    double lowLimit = 0.92 * rrAverage2;
    double highLimit = 1.16 * rrAverage2;

    //sudden rate change: the newest interval is outside the limits,
    //so RR AVERAGE1 is used for the missed-beat limit
    bool suddenChange = rrInterval < lowLimit || rrInterval > highLimit;

    //regularity can only be judged once there are 8 recent intervals
    if (recentRR.size() == 8) {
        bool regular = true;
        for (int value : recentRR) {
            if (value < lowLimit || value > highLimit) {
                regular = false;
            }
        }
        irregularRhythm = !regular;

        //for a regular rhythm the paper sets RR AVERAGE2 = RR AVERAGE1,
        //so the regular interval set becomes the recent interval set
        if (regular) {
            regularRR = recentRR;
            rrAverage2 = rrAverage1;
        }
    }

    if (suddenChange) {
        rrMissedLimit = 1.66 * rrAverage1;
    } else {
        rrMissedLimit = 1.66 * rrAverage2;
    }
}



/* checks if a peak that comes soon after a QRS is a T-wave
   if the peak is less than 360 ms after the last QRS and its slope is
   less than half of the slope of that QRS, it is a T-wave. */
bool PanTompkinsDetector::isTWave(const PeakCandidate& candidate) {
    if (lastQRSPosition == -1) {
        return false;
    }

    int tWaveWindow = static_cast<int>(0.36 * processingRate);   //360 ms
    bool tooClose = candidate.position - lastQRSPosition < tWaveWindow;
    bool smallSlope = candidate.slope < 0.5 * lastQRSSlope;

    return tooClose && smallSlope;
}



//saves a peak as a QRS complex and updates the RR information. 
void PanTompkinsDetector::acceptQRS(const PeakCandidate& candidate) {
    qrsPositions.push_back(candidate.positionF);

    if (lastQRSPosition != -1) {
        updateRRAverages(candidate.position - lastQRSPosition);
    }

    lastQRSPosition = candidate.position;
    lastQRSSlope = candidate.slope;

    //noise peaks before this QRS are no longer needed for search-back
    vector<PeakCandidate> laterPeaks;
    for (const PeakCandidate& peak : noisePeaks) {
        if (peak.position > candidate.position) {
            laterPeaks.push_back(peak);
        }
    }
    noisePeaks = laterPeaks;
}



//updates the noise levels and keeps the peak for a possible search-back
void PanTompkinsDetector::markAsNoise(const PeakCandidate& candidate) {
    noisePeakI = 0.125 * candidate.peakI + 0.875 * noisePeakI;
    noisePeakF = 0.125 * candidate.peakF + 0.875 * noisePeakF;

    noisePeaks.push_back(candidate);
}



/* search-back for a missed beat
   from the noise peaks since the last QRS, takes the highest one that is
   above THRESHOLD2 in both signals and treats it as a QRS
   ->> return true if a missed beat was found */
bool PanTompkinsDetector::searchBack() {
    int refractoryPeriod = static_cast<int>(0.20 * processingRate);
    int bestIndex = -1;

    for (int k = 0; k < noisePeaks.size(); k++) {
        const PeakCandidate& peak = noisePeaks[k];

        bool afterRefractory = peak.position - lastQRSPosition > refractoryPeriod;
        bool aboveThresholdI = peak.peakI > thresholdI2;
        bool aboveThresholdF = peak.peakF > thresholdF2;

        if (afterRefractory && aboveThresholdI && aboveThresholdF && !isTWave(peak)) {
            if (bestIndex == -1 || peak.peakI > noisePeaks[bestIndex].peakI) {
                bestIndex = k;
            }
        }
    }

    if (bestIndex == -1) {
        return false;
    }

    PeakCandidate missedBeat = noisePeaks[bestIndex];

    //search-back uses a different update from the paper (0.25 / 0.75)
    signalPeakI = 0.25 * missedBeat.peakI + 0.75 * signalPeakI;
    signalPeakF = 0.25 * missedBeat.peakF + 0.75 * signalPeakF;

    acceptQRS(missedBeat);
    updateThresholds();

    return true;
}




/* finds the exact R-peak near the detected QRS in the ECG signal
   (the paper lets you choose the fiducial point, this method is my choice)
   - the local mean is removed so the baseline does not affect the result
   return R-peak position */
int PanTompkinsDetector::findRPeak(const vector<double>& signal, int approximatePosition) {
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