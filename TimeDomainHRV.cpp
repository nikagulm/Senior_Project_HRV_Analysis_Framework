#include "TimeDomainHRV.h"
#include <cmath>
using namespace std;

//nn intervals in seconds & results in milliseconds
TimeDomainResult TimeDomainHRV::calculate(const vector<double>& nn) {
    TimeDomainResult result{};
    int n = nn.size();
    if (n < 2) {
        return result;
    }

    //mean NN
    double sum = 0.0;
    for (double x : nn) {
        sum += x;
    }
    double mean = sum / n;

    //standard deviation of the NN intervals // SSDN
    double squares = 0.0;
    for (double x : nn) {
        squares += (x - mean) * (x - mean);
    }
    double sdnn = sqrt(squares / (n - 1));

    //RMSSD and pNN50 from the differences between successive intervals
    double diffSquares = 0.0;
    int over50 = 0;
    for (int i = 1; i < n; i++) {
        double diff = nn[i] - nn[i - 1];
        diffSquares += diff * diff;
        if (abs(diff) > 0.050) {
            over50++;
        }
    }

    result.meanNN = mean * 1000.0;
    result.sdnn = sdnn * 1000.0;
    result.rmssd = sqrt(diffSquares / (n - 1)) * 1000.0;
    result.pnn50 = 100.0 * over50 / (n - 1);

    return result;
}