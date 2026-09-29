using namespace std;

#ifndef PANTOMPKINSDETECTOR_H
#define PANTOMPKINSDETECTOR_H

#include <vector>

class PanTompkinsDetector {
public:
    vector<int> detect(const vector<double>& ecgSignal);

private:
    vector<double> bandPassFilter(const vector<double>& signal);
    vector<double> derivative(const vector<double>& signal);
    vector<double> squareSignal(const vector<double>& signal);
    vector<double> movingWindowIntegration(const vector<double>& signal);
    vector<int> detectPeaks(const vector<double>& signal);
};

#endif