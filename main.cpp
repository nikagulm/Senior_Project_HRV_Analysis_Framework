#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include "headerParser.h"
#include "dataParser.h"
#include "annotationParser.h"
#include "detectorEvaluator.h"
#include "RPeakDetector.h"
#include "PanTompkinsDetector.h"
#include "HamiltonDetector.h"
#include "HilbertDetector.h"

using namespace std;

int main() {
    vector<string> records = {"100", "101"};

    for (const string& record : records) {
        //read the record
        HeaderParser header;
        header.parse(record + ".hea");
        DataParser data;
        data.parseFormat212(record + ".dat");
        vector<int> channel1 = data.getChannel1();
        AnnotationParser annotations;
        vector<int> reference = annotations.parse(record + ".atr");

        //ADC values -> millivolts
        const SignalInfo& info = header.getSignals()[0];
        vector<double> ecg(channel1.size());
        for (int i = 0; i < channel1.size(); i++) {
            ecg[i] = (channel1[i] - info.adcZero) / static_cast<double>(info.gain);
        }

        int fs = header.getSamplingRate();
        PanTompkinsDetector panTompkins(fs);
        HamiltonDetector hamilton(fs);
        HilbertDetector hilbert(fs);

        //all detectors through the same interface
        vector<pair<string, RPeakDetector*>> detectors = {
            {"Pan-Tompkins", &panTompkins},
            {"Hamilton", &hamilton},
            {"Hilbert", &hilbert}
        };

        DetectorEvaluator evaluator;
        int tolerance = static_cast<int>(round(0.150 * fs)); //150 ms

        cout << "\nRecord " << record << " (" << reference.size() << " reference beats)" << endl;
        for (auto& d : detectors) {
            vector<int> peaks = d.second->detect(ecg);
            EvaluationResult r = evaluator.evaluate(peaks, reference, tolerance);
            cout << "  " << d.first << ": detected=" << peaks.size()
                 << " TP=" << r.truePositives << " FP=" << r.falsePositives << " FN=" << r.falseNegatives
                 << " Se=" << r.sensitivity * 100 << "% PPV=" << r.ppv * 100 << "%" << endl;
        }
    }
    return 0;
}