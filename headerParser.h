#ifndef HEADERPARSER_H
#define HEADERPARSER_H
#include <string>
#include <vector>

struct SignalInfo {
    std::string dataFileName;
    int format;
    int gain;
    int adcResolution;
    int adcZero;
    int initialValue;
    int checksum;
    int blockSize;
    std::string signalName;
};

class HeaderParser {
public:
    void parse(const std::string& fileName);
    std::string getRecordName() const;
    int getNumberOfSignals() const;
    int getSamplingRate() const;
    int getNumberOfSamples() const;
    const std::vector<SignalInfo>& getSignals() const;

private:
    std::string recordName;
    int numberOfSignals;
    int samplingRate;
    int numberOfSamples;

    std::vector<SignalInfo> signals;
};


#endif