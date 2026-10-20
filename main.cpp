#include <iostream>
#include <vector>
#include <string>

#include "PanTompkinsDetector.h"
#include "HamiltonDetector.h"
#include "headerParser.h"
#include "dataParser.h"
#include "annotationParser.h"
#include "detectorEvaluator.h"

using namespace std;

int main() {
    string record = "100";

    //parse the header file
    HeaderParser parser;
    parser.parse(record + ".hea");

    cout << "Sampling rate: " << parser.getSamplingRate() << endl;
    cout << "Stored signals: " << parser.getSignals().size() << endl;

    //parse ECG samples
    DataParser dataParser;
    vector<int> samples = dataParser.parseFormat212(record + ".dat");

    cout << "Decoded samples: " << samples.size() << endl;
    cout << "Channel 1 samples: " << dataParser.getChannel1().size() << endl;
    cout << "Channel 2 samples: " << dataParser.getChannel2().size() << endl;

    //use the first ECG channel
    vector<int> channel1 = dataParser.getChannel1();
    vector<double> ecgSignal(channel1.begin(), channel1.end());

    //parse reference heartbeat annotations
    AnnotationParser annotationParser;
    vector<int> annotations = annotationParser.parse(record + ".atr");

    cout << "Annotations: " << annotations.size() << endl;


    //run Pan-Tompkins
    PanTompkinsDetector panTompkins(parser.getSamplingRate());
    vector<int> panTompkinsPeaks = panTompkins.detect(ecgSignal);

    cout << "\nPan-Tompkins detected peaks: " << panTompkinsPeaks.size() << endl;

    //compare the first Pan-Tompkins peaks with the reference annotations
    cout << "\nFirst Pan-Tompkins peaks vs annotations:" << endl;

    for (int i = 0; i < 10 && i < panTompkinsPeaks.size() && i < annotations.size(); i++) {
        cout << "Detected: " << panTompkinsPeaks[i] << "   Reference: " << annotations[i] << endl;
        cout << "Difference: " << panTompkinsPeaks[i] - annotations[i] << endl;
    }


    //run Hamilton
    HamiltonDetector hamilton(parser.getSamplingRate());
    vector<int> hamiltonPeaks = hamilton.detect(ecgSignal);

    cout << "\nHamilton detected peaks: " << hamiltonPeaks.size() << endl;

    //compare the first Hamilton peaks with the reference annotations
    cout << "\nFirst Hamilton peaks vs annotations:" << endl;

    for (int i = 0; i < 10 && i < hamiltonPeaks.size() && i < annotations.size(); i++) {
        cout << "Detected: " << hamiltonPeaks[i] << "   Reference: " << annotations[i] << endl;
        cout << "Difference: " << hamiltonPeaks[i] - annotations[i] << endl;
    }


    //evaluate both detectors using the same tolerance windows
    DetectorEvaluator evaluator;
    vector<int> tolerances = {15, 20, 25, 30};

    cout << "\n--- Pan-Tompkins Evaluation ---" << endl;

    for (int tolerance : tolerances) {
        EvaluationResult result = evaluator.evaluate(panTompkinsPeaks, annotations, tolerance);

        cout << "\nTolerance: " << tolerance << " samples" << endl;
        cout << "True positives: " << result.truePositives << endl;
        cout << "False positives: " << result.falsePositives << endl;
        cout << "False negatives: " << result.falseNegatives << endl;
        cout << "Sensitivity: " << result.sensitivity << endl;
        cout << "PPV: " << result.ppv << endl;
    }

    cout << "\n--- Hamilton Evaluation ---" << endl;
    for (int tolerance : tolerances) {
        EvaluationResult result = evaluator.evaluate(hamiltonPeaks, annotations, tolerance);

        cout << "\nTolerance: " << tolerance << " samples" << endl;
        cout << "True positives: " << result.truePositives << endl;
        cout << "False positives: " << result.falsePositives << endl;
        cout << "False negatives: " << result.falseNegatives << endl;
        cout << "Sensitivity: " << result.sensitivity << endl;
        cout << "PPV: " << result.ppv << endl;
    }

    return 0;
}