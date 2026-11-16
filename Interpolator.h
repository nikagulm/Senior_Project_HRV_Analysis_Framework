#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include <vector>

//turns the uneven NN series into an evenly sampled signal for the spectrum
class Interpolator {
public:
    static std::vector<double> cubicSpline(const std::vector<double>& times, const std::vector<double>& values, double fs);
    static void detrend(std::vector<double>& signal);
};

#endif