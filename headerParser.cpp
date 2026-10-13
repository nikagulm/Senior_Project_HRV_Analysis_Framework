#include "headerParser.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

void HeaderParser::parse(const string& fileName) {
    ifstream file(fileName);
    if (!file.is_open()) {
        cout << "Can't open the file" << endl;
        return;
    }

    string firstLine;
    getline(file, firstLine);

    istringstream lineStream(firstLine);
    lineStream >> recordName
               >> numberOfSignals
               >> samplingRate
               >> numberOfSamples;

    cout << "Record name: " << recordName << endl;
    cout << "Number of signals: " << numberOfSignals << endl;
    cout << "Sampling rate: " << samplingRate << endl;
    cout << "Number of samples: " << numberOfSamples << endl;

    string signalLine;
    for (int i = 0; i < numberOfSignals; i++) {
        getline(file, signalLine);
        istringstream signalStream(signalLine);

        string dataFileName;
        int format, gain, adcResolution, adcZero, initialValue, checksum, blockSize;
        string signalName;

        //extract value from each signal line
        signalStream >> dataFileName
                    >> format
                    >> gain
                    >> adcResolution
                    >> adcZero
                    >> initialValue
                    >> checksum
                    >> blockSize
                    >> signalName;
        
        SignalInfo info;
        info.dataFileName = dataFileName;
        info.format = format;
        info.gain = gain;
        info.adcResolution = adcResolution;
        info.adcZero = adcZero;
        info.initialValue = initialValue;
        info.checksum = checksum;
        info.blockSize = blockSize;
        info.signalName = signalName;

        signals.push_back(info);
        cout << "Data file: " << dataFileName << endl;
        cout << "Format: " << format << endl;
        cout << "Gain: " << gain << endl;
        cout << "ADC resolution: " << adcResolution << endl;
        cout << "ADC zero: " << adcZero << endl;
        cout << "Initial value: " << initialValue << endl;
        cout << "Checksum: " << checksum << endl;
        cout << "Block size: " << blockSize << endl;
        cout << "Signal name: " << signalName << endl;    
        cout << "Signal line " << i + 1 << ": " << signalLine << endl;
    }
    cout << "Stored signals: " << signals.size() << endl;

    //file.close();
}

string HeaderParser::getRecordName() const {
    return recordName;
}

int HeaderParser::getNumberOfSignals() const {
    return numberOfSignals;
}

int HeaderParser::getSamplingRate() const {
    return samplingRate;
}

int HeaderParser::getNumberOfSamples() const {
    return numberOfSamples;
}

const vector<SignalInfo>& HeaderParser::getSignals() const {
    return signals;
}