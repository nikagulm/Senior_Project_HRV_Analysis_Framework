#include "RRIntervalProcessor.h"
using namespace std;

RRIntervalProcessor::RRIntervalProcessor(int samplingRate) {
    this->samplingRate = samplingRate;
}

vector<double> RRIntervalProcessor::rrIntervals(const vector<int>& rPeaks) { //convert R-peak positions into RR intervals 
    vector<double> rr;
    for (int i = 1; i < rPeaks.size(); i++) {
        rr.push_back(static_cast<double>(rPeaks[i] - rPeaks[i - 1]) / samplingRate);
    }
    return rr;
}

//remove intervals caused by false beats 
vector<double> RRIntervalProcessor::nnIntervals(const vector<double>& rr) {
    vector<double> nn;
    for (int i = 0; i < rr.size(); i++) {
        if (rr[i] < minInterval || rr[i] > maxInterval) { //impossible values
            continue;
        }

        if (i > 0 && abs(rr[i] - rr[i - 1]) > maxChange * rr[i - 1]) {
            continue;
        }
        nn.push_back(rr[i]);
        
    }
    return nn; //result are clean NN intervals
}