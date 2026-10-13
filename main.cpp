#include <iostream>
#include <vector>

#include "PanTompkinsDetector.h"
#include "headerParser.h"
#include "dataParser.h"
#include "annotationParser.h"

using namespace std;

int main() {
    // Parse the header file
    HeaderParser parser;
    parser.parse("100.hea");

    cout << "Testing getter sampling rate: "
         << parser.getSamplingRate() << endl;

    cout << "Testing getter stored signals: "
         << parser.getSignals().size() << endl;


    // Parse the ECG data file
    DataParser dataParser;
    vector<int> samples = dataParser.parseFormat212("100.dat");

    cout << "Decoded samples: " << samples.size() << endl;

    for (int i = 0; i < 10 && i < samples.size(); i++) {
        cout << samples[i] << " ";
    }

    cout << endl;


    // Get the two ECG channels
    vector<int> channel1 = dataParser.getChannel1();

    cout << "Channel 1 samples: "
         << dataParser.getChannel1().size() << endl;

    cout << "Channel 2 samples: "
         << dataParser.getChannel2().size() << endl;


    // Convert channel 1 to double values for Pan-Tompkins
    vector<double> ecgSignal(channel1.begin(), channel1.end());

    // Run Pan-Tompkins on the real ECG signal
    PanTompkinsDetector detector(parser.getSamplingRate());

    vector<int> detectedPeaks = detector.detect(ecgSignal);

    cout << "Detected peaks: "
         << detectedPeaks.size() << endl;


    // Parse the reference annotations
    AnnotationParser annotationParser;
    vector<int> annotations = annotationParser.parse("100.atr");

    cout << "Annotations: "
         << annotations.size() << endl;

    for (int i = 0; i < 10 && i < annotations.size(); i++) {
        cout << annotations[i] << " ";
    }

    cout << "First detected peaks vs annotations:" << endl;

    for (int i = 0; i < 10 && i < detectedPeaks.size() && i < annotations.size(); i++) {
        cout << "Detected: " << detectedPeaks[i] << "   Reference: " << annotations[i] << endl;
        cout << "Difference: " << detectedPeaks[i] - annotations[i] << endl;
    }
    cout << endl;

    return 0;
}