#ifndef TIMEDOMAINHRV_H
#define TIMEDOMAINHRV_H
#include <vector>

struct TimeDomainResult {
    double meanNN;
    double sdnn;
    double rmssd;
    double pnn50;
};

class TimeDomainHRV {
public:
    TimeDomainResult calculate(const std::vector<double>& nn);
};

#endif