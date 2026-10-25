#include "HRVAnalysisPipeline.h"
#include <iostream>
#include <cmath>

using namespace std;

void HRVAnalysisPipeline::run(const vector<double>& ecgSignal, const vector<int>& annotations, int samplingRate, const vector<pair<string, RPeakDetector*>>& detectors) {
    DetectorEvaluator evaluator;
    int standardTolerance = static_cast<int>(round(0.150 * samplingRate));
    vector<int> tolerances = {15, 20, 25, 30, standardTolerance};

    for (const auto& entry : detectors) {
        const string& detectorName = entry.first;
        RPeakDetector* detector = entry.second;

        vector<int> detectedPeaks = detector->detect(ecgSignal);

        cout << "\n--- " << detectorName << " ---" << endl;
        cout << "Detected peaks: " << detectedPeaks.size() << endl;
        cout << "\nFirst detected peaks vs annotations:" << endl;
        for (size_t i = 0; i < 10 && i < detectedPeaks.size() && i < annotations.size(); i++) {
            cout << "Detected: " << detectedPeaks[i] << "   Reference: " << annotations[i] << endl;
            cout << "Difference: " << detectedPeaks[i] - annotations[i] << endl;
        }

        cout << "\n--- " << detectorName << " Evaluation ---" << endl;
        for (int tolerance : tolerances) {
            EvaluationResult result = evaluator.evaluate(detectedPeaks, annotations, tolerance);

            cout << "\nTolerance: " << tolerance << " samples" << endl;
            cout << "True positives: " << result.truePositives << endl;
            cout << "False positives: " << result.falsePositives << endl;
            cout << "False negatives: " << result.falseNegatives << endl;
            cout << "Sensitivity: " << result.sensitivity << endl;
            cout << "PPV: " << result.ppv << endl;
        }
    }
}
