#include <iostream>
#include <vector>
#include <string>
#include "PanTompkinsDetector.h"
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

    //run Pan-Tompkins
    PanTompkinsDetector detector(parser.getSamplingRate());
    vector<int> detectedPeaks = detector.detect(ecgSignal);

    cout << "Detected peaks: " << detectedPeaks.size() << endl;

    //parse reference heartbeat annotations
    AnnotationParser annotationParser;
    vector<int> annotations = annotationParser.parse(record + ".atr");

    cout << "Annotations: " << annotations.size() << endl;

    //compare the first detected peaks with the reference annotations
    cout << "\nFirst detected peaks vs annotations:" << endl;

    for (int i = 0; i < 10 && i < detectedPeaks.size() && i < annotations.size(); i++) {
        cout << "Detected: " << detectedPeaks[i] << "   Reference: " << annotations[i] << endl;
        cout << "Difference: " << detectedPeaks[i] - annotations[i] << endl;
    }

    //evaluate detector performance using different tolerance windows
    DetectorEvaluator evaluator;
    vector<int> tolerances = {15, 20, 25, 30};

    for (int tolerance : tolerances) {
        EvaluationResult result = evaluator.evaluate(detectedPeaks, annotations, tolerance);

        cout << "\nTolerance: " << tolerance << " samples" << endl;
        cout << "True positives: " << result.truePositives << endl;
        cout << "False positives: " << result.falsePositives << endl;
        cout << "False negatives: " << result.falseNegatives << endl;
        cout << "Sensitivity: " << result.sensitivity << endl;
        cout << "PPV: " << result.ppv << endl;
    }

    return 0;
}