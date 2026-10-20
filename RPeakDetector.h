#ifndef RPEAKDETECTOR_H
#define RPEAKDETECTOR_H
#include <vector>

class RPeakDetector {
public:
    virtual ~RPeakDetector() = default;
    virtual std::vector<int> detect(const std::vector<double>& ecgSignal) = 0;
    };

#endif
