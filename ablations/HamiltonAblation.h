#ifndef HAMILTONABLATION_H
#define HAMILTONABLATION_H
#include "../HamiltonDetector.h"
#include "AblationVariant.h"

//Hamilton with exactly one stage replaced, everything else is the reference code
class HamiltonAblation : public HamiltonDetector {
public:
    HamiltonAblation(int samplingRate, AblationVariant variant);

    std::vector<int> detect(const std::vector<double>& ecgSignal) override;

protected:
    AblationVariant variant;
    bool thresholdFrozen = false;
    int averageLength = 9; //45 ms at 200 Hz

    std::vector<double> lowPassFilter(const std::vector<double>& signal) override;
    std::vector<double> highPassFilter(const std::vector<double>& signal) override;
    void updateThreshold() override;
    bool useSearchBack() const override;
};

#endif