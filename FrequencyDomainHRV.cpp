#include "FrequencyDomainHRV.h"
#include "Interpolator.h"
using namespace std;

FrequencyDomainHRV::FrequencyDomainHRV(SpectrumEstimator* estimator) {
    this->estimator = estimator; 
}

FrequencyDomainResult FrequencyDomainHRV::calculate(const vector<double>& nnTimes, const vector<double>& nn) {
    FrequencyDomainResult result = {0.0, 0.0, 0.0, 0.0, 0};

    vector<double> nnMs;
    for (double x : nn) {
        nnMs.push_back(x * 1000.0); //NN in ms
    }

    //even samples at 4 Hz
    vector<double> signal = Interpolator::cubicSpline(nnTimes, nnMs, resampleRate);
    int n = segmentLength;
    int m = 1;
    while (m < n) {
        m *= 2;
    }
    double binWidth = resampleRate / m;

    //all frequencies from VLF up to the end of HF
    vector<double> frequencies;
    for (int k = 1; k * binWidth <= hfHigh; k++) {
        if (k * binWidth >= vlfLow) {
            frequencies.push_back(k * binWidth);
        }
    }

    //5 minute segments w band powers averaged over them
    for (int start = 0; start + n <= signal.size(); start += n) {
        vector<double> segment(signal.begin() + start, signal.begin() + start + n);
        Interpolator::detrend(segment); //remove mean and drift

        vector<double> power = estimator->power(segment, resampleRate, frequencies);

        //add each frequency power to its band
        double vlf = 0.0, lf = 0.0, hf = 0.0;
        for (int i = 0; i < frequencies.size(); i++) {
            double p = 2.0 * power[i] / (static_cast<double>(n) * m); //convert to ms^2
            double f = frequencies[i];
            if (f < lfLow) {
                vlf += p;
            } else if (f < hfLow) {
                lf += p;
            } else {
                hf += p;
            }
        }
        result.vlf += vlf;
        result.lf += lf;
        result.hf += hf;
        result.segments++;
    }

    //average over segments
    if (result.segments > 0) {
        result.vlf /= result.segments;
        result.lf /= result.segments;
        result.hf /= result.segments;
        result.lfhf = result.lf / result.hf;
    }
    return result;
}