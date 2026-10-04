#include <iostream>
#include <vector>
#include "PanTompkinsDetector.h"
using namespace std;

int main() {
    int samplingRate = 200;
    PanTompkinsDetector detector(samplingRate);
    vector<double> signal(1000, 0.0);

    signal[100] = 10.0;
    signal[300] = 10.0;
    signal[500] = 10.0;
    signal[700] = 10.0;
    signal[900] = 10.0;

    vector<int> peaks = detector.detect(signal);
    cout << "Peaks found are: " << endl;
    for (int peak : peaks) {
        cout << peak << endl;
    }

    return 0;
}