#ifndef PANTOMPKINSABLATION_H
#define PANTOMPKINSABLATION_H

#include "../PanTompkinsDetector.h"
#include "AblationVariant.h"

//Pan-Tompkins with exactly one stage replaced, everything else is the reference code
class PanTompkinsAblation : public PanTompkinsDetector {
public:
    PanTompkinsAblation(int samplingRate, AblationVariant variant);

    std::vector<int> detect(const std::vector<double>& ecgSignal) override;

protected:
    AblationVariant variant;
    bool thresholdsFrozen = false;
    int averageLength = 9; //45 ms at 200 Hz

    std::vector<double> lowPassFilter(const std::vector<double>& signal) override;
    std::vector<double> highPassFilter(const std::vector<double>& signal) override;
    int learnThresholds(const std::vector<double>& integrated, const std::vector<double>& bandPassed) override;
    void updateThresholds() override;
    bool searchBack() override;
};

#endif