#include "RRIntervalProcessor.h"
#include <cmath>
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


//checks if one interval is normal
bool RRIntervalProcessor::isNormal(const vector<double>& rr, int i) {
    if (rr[i] < minInterval || rr[i] > maxInterval) { //impossible values
        return false;
    }
    if (i > 0 && abs(rr[i] - rr[i - 1]) > maxChange * rr[i - 1]) { //sudden jump
        return false;
    }
    return true; 
}

//remove intervals caused by false beats 
vector<double> RRIntervalProcessor::nnIntervals(const vector<double>& rr) {
    vector<double> nn;
    for (int i = 0; i < rr.size(); i++) {
        if (isNormal(rr, i)) {
            nn.push_back(rr[i]);
        }
    }
    return nn;


    /* for (int i = 0; i < rr.size(); i++) {
        if (rr[i] < minInterval || rr[i] > maxInterval) { //impossible values
            continue;
        }

        if (i > 0 && abs(rr[i] - rr[i - 1]) > maxChange * rr[i - 1]) {
            continue;
        }
        nn.push_back(rr[i]);
    
    }
    return nn; //result are clean NN intervals */

}

//time of each NN beat, needed for the spline
vector<double> RRIntervalProcessor::nnTimes(const vector<double>& rr) {
    vector<double> times;
    double t = 0.0;
    for (int i = 0; i < rr.size(); i++) {
        t += rr[i]; //count every interval so the gaps stay in the right place
        if (isNormal(rr, i)) {
            times.push_back(t);
        }
    }
    return times;
}