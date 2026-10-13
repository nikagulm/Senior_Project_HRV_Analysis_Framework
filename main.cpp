#include <iostream>
#include <vector>
#include "PanTompkinsDetector.h"
#include "headerParser.h"
#include "dataParser.h"
#include "annotationParser.h"
#include "detectorEvaluator.h"

using namespace std;

int main() {
    //parse the header file
    HeaderParser parser;
    parser.parse("100.hea");

    cout << "Testing getter sampling rate: "
         << parser.getSamplingRate() << endl;

    cout << "Testing getter stored signals: "
         << parser.getSignals().size() << endl;


    //parse the ECG data file
    DataParser dataParser;
    vector<int> samples = dataParser.parseFormat212("100.dat");

    cout << "Decoded samples: " << samples.size() << endl;

    for (int i = 0; i < 10 && i < samples.size(); i++) {
        cout << samples[i] << " ";
    }

    cout << endl;


    //get channel 1
    vector<int> channel1 = dataParser.getChannel1();

    cout << "Channel 1 samples: "
         << dataParser.getChannel1().size() << endl;

    cout << "Channel 2 samples: "
         << dataParser.getChannel2().size() << endl;


    //convert channel 1 to double values for Pan-Tompkins
    vector<double> ecgSignal(channel1.begin(), channel1.end());


    //run Pan-Tompkins
    PanTompkinsDetector detector(parser.getSamplingRate());

    vector<int> detectedPeaks = detector.detect(ecgSignal);

    cout << "Detected peaks: "
         << detectedPeaks.size() << endl;

    vector<int> correctedPeaks = detectedPeaks;
    for (int& peak : correctedPeaks) {
        peak -= 29;
    }


    //parse the reference annotations
    AnnotationParser annotationParser;
    vector<int> annotations = annotationParser.parse("100.atr");

    cout << "Annotations: "
         << annotations.size() << endl;

    for (int i = 0; i < 10 && i < annotations.size(); i++) {
        cout << annotations[i] << " ";
    }

    cout << endl;

    //compare the first detected peaks with reference annotations
    cout << "\nFirst detected peaks vs annotations:" << endl;
    for (int i = 0;
         i < 10 &&
         i < detectedPeaks.size() &&
         i < annotations.size();
         i++) {

        cout << "Detected: " << detectedPeaks[i]
             << "   Reference: " << annotations[i] << endl;

        cout << "Difference: "
             << detectedPeaks[i] - annotations[i]
             << endl;
    }


    //evaluate the detector using different tolerance windows
    DetectorEvaluator evaluator;

    vector<int> tolerances = {15, 20, 25, 30};
    for (int tolerance : tolerances) {
        EvaluationResult result = evaluator.evaluate(
            correctedPeaks,
            annotations,
            tolerance
        );

        cout << "\nTolerance: " << tolerance << " samples" << endl;
        cout << "True positives: " << result.truePositives << endl;
        cout << "False positives: " << result.falsePositives << endl;
        cout << "False negatives: " << result.falseNegatives << endl;
        cout << "Sensitivity: " << result.sensitivity << endl;
        cout << "PPV: " << result.ppv << endl;
    }

    return 0;
}