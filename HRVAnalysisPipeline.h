#ifndef HRVANALYSISPIPELINE_H
#define HRVANALYSISPIPELINE_H
#include <vector>
#include <string>
#include <utility>
#include "RPeakDetector.h"
#include "detectorEvaluator.h"

class HRVAnalysisPipeline {
public:
    void run(const std::vector<double>& ecgSignal, const std::vector<int>& annotations, int samplingRate, const std::vector<std::pair<std::string, RPeakDetector*>>& detectors);
};

#endif
