#include "Interpolator.h"
using namespace std;

//natural cubic spline through the beats read every 1/fs seconds
vector<double> Interpolator::cubicSpline(const vector<double>& times, const vector<double>& values, double fs) {
    int n = times.size();
    vector<double> result;
    if (n < 3) {
        return result;
    }

    //time between beats
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; i++) {
        h[i] = times[i + 1] - times[i];
    }

    vector<double> M(n, 0.0); //bend of the curve at each beat 
    vector<double> diag(n, 0.0), rhs(n, 0.0); //one equation per beat so the curve has no corners
    for (int i = 1; i < n - 1; i++) {
        diag[i] = 2.0 * (h[i - 1] + h[i]);
        rhs[i] = 6.0 * ((values[i + 1] - values[i]) / h[i] - (values[i] - values[i - 1]) / h[i - 1]);
    }

    //Thomas alg
    for (int i = 2; i < n - 1; i++) {
        double w = h[i - 1] / diag[i - 1];
        diag[i] -= w * h[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }

    //back substitution
    for (int i = n - 2; i >= 1; i--) {
        M[i] = (rhs[i] - h[i] * M[i + 1]) / diag[i];
    }

    //read the curve at equal time steps
    int j = 0;
    for (double t = times[0]; t <= times[n - 1]; t += 1.0 / fs) {
        //find the two beats around t
        while (j < n - 2 && t > times[j + 1]) {
            j++;
        }
        double a = times[j + 1] - t;
        double b = t - times[j];
        double hj = h[j];

        //spline value between beat j and j+1
        double value = M[j] * a * a * a / (6.0 * hj) + M[j + 1] * b * b * b / (6.0 * hj) + (values[j] / hj - M[j] * hj / 6.0) * a + (values[j + 1] / hj - M[j + 1] * hj / 6.0) * b;
        result.push_back(value);
    }
    return result;
}

//removes the mean and the linear trend before the spectrum
void Interpolator::detrend(vector<double>& signal) {
    int n = signal.size();
    if (n < 2) {
        return;
    }

    double meanI = (n - 1) / 2.0;
    double meanX = 0.0;
    for (double x : signal) {
        meanX += x;
    }
    meanX /= n;

    //slope of the best-fit line
    double top = 0.0, bottom = 0.0;
    for (int i = 0; i < n; i++) {
        top += (i - meanI) * (signal[i] - meanX);
        bottom += (i - meanI) * (i - meanI);
    }
    double slope = top / bottom;

    //subtract the line
    for (int i = 0; i < n; i++) {
        signal[i] -= meanX + slope * (i - meanI);
    }
}