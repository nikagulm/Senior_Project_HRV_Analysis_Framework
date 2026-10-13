#ifndef DETECTOREVALUATOR_H
#define DETECTOREVALUATOR_H

#include <vector>

struct EvaluationResult {
    int truePositives;
    int falsePositives;
    int falseNegatives;
    double sensitivity;
    double ppv;
};

class DetectorEvaluator {
public:
    EvaluationResult evaluate(
        const std::vector<int>& detectedPeaks,
        const std::vector<int>& annotations,
        int toleranceSamples
    );
};

#endif