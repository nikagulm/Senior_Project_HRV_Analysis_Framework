//time-domain HRV, RR/NN intervals, spline and detrending
#include "TestUtils.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include "Interpolator.h"
#include "RRIntervalProcessor.h"
#include "TimeDomainHRV.h"
using namespace std;

void testTimeDomain() {
    printf("\nTime-domain HRV\n");
    //values calculated by hand differences +60, -70, +20, +90 ms
    TimeDomainHRV td;
    TimeDomainResult r = td.calculate({0.80, 0.86, 0.79, 0.81, 0.90});
    check("mean NN (expected 832.0 ms)", abs(r.meanNN - 832.0) < 0.01, fmt("%.2f ms", r.meanNN));
    check("SDNN (expected 46.58 ms)", abs(r.sdnn - 46.58) < 0.01, fmt("%.2f ms", r.sdnn));
    check("RMSSD (expected 65.19 ms)", abs(r.rmssd - 65.19) < 0.01, fmt("%.2f ms", r.rmssd));
    check("pNN50 (expected 75.00 %)", abs(r.pnn50 - 75.0) < 0.01, fmt("%.2f %%", r.pnn50));
    TimeDomainResult empty = td.calculate({0.8});
    check("a single interval gives zeros instead of a division by zero", empty.sdnn == 0 && empty.rmssd == 0, "");
}

void testNNRules() {
    printf("\nRR and NN intervals\n");
    RRIntervalProcessor rrp(360);
    vector<double> rr = rrp.rrIntervals({0, 288, 576, 900});
    check("R-peaks 288 samples apart at 360 Hz give 0.8 s", rr.size() == 3 && abs(rr[0] - 0.8) < 1e-12 && abs(rr[2] - 0.9) < 1e-12,
          fmt("%.3f, %.3f, %.3f s", rr[0], rr[1], rr[2]));

    //0.25s and 2.5s are impossible,
    //1.10s is a jump of more than 20 %
    vector<double> raw = {0.80, 0.82, 0.25, 0.80, 0.84, 1.10, 0.86, 2.50, 0.85};
    vector<double> nn = rrp.nnIntervals(raw);
    string kept;
    for (double x : nn) {
        kept += fmt("%.2f ", x);
    }
    //kept 0.80, 0.82 and 0.84 (0.84 differs from 0.80 by 5 %) 0.86 is rejected after the 1.10 jump
    check("range 0.3-2.0 s and 20 % rule (expected 0.80 0.82 0.84)", nn.size() == 3 && abs(nn[2] - 0.84) < 1e-12, "kept " + kept);
    vector<double> times = rrp.nnTimes(raw);
    check("NN times keep the gaps of rejected intervals (expected 0.80 1.62 3.51 s)",
          times.size() == 3 && abs(times[2] - 3.51) < 1e-9, fmt("%.2f, %.2f, %.2f s", times[0], times[1], times[2]));
}

void testSpline() {
    printf("\nSpline and detrending\n");
    //a natural cubic spline reproduces a straight line exactly even with uneven knots
    vector<double> t = {0.0, 0.7, 1.5, 2.6, 3.2, 4.0};
    vector<double> y;
    for (double x : t) {
        y.push_back(2.0 * x + 5.0);
    }
    vector<double> s = Interpolator::cubicSpline(t, y, 4.0);
    double err = 0;
    for (int i = 0; i < s.size(); i++) {
        err = max(err, abs(s[i] - (2.0 * i / 4.0 + 5.0)));
    }
    check("spline through a straight line, read at 4 Hz", err < 1e-9 && s.size() == 17, fmt("%.0f samples, max error %.1e", s.size(), err));

    //the spline must pass through every knot of a curved series
    vector<double> t2 = {0.0, 0.8, 1.6, 2.4, 3.2, 4.0};
    vector<double> y2 = {800, 850, 790, 820, 880, 810};
    vector<double> s2 = Interpolator::cubicSpline(t2, y2, 4.0); //knots fall on every 0.8 s = 3.2 samples
    double knotErr = max({abs(s2[0] - 800), abs(s2[16] - 810)});
    check("spline passes through the first and last knot", knotErr < 1e-9, fmt("max error %.1e", knotErr));
    vector<double> trend(1200);
    for (int i = 0; i < 1200; i++) {
        trend[i] = 3.0 + 0.5 * i + 10 * sin(2 * M_PI * i / 40.0);
    }
    Interpolator::detrend(trend);
    double mean = 0, slopeTop = 0, slopeBottom = 0;
    for (int i = 0; i < 1200; i++) {
        mean += trend[i] / 1200;
    }
    for (int i = 0; i < 1200; i++) {
        slopeTop += (i - 599.5) * trend[i];
        slopeBottom += (i - 599.5) * (i - 599.5);
    }
    check("detrending removes mean and linear trend", abs(mean) < 1e-9 && abs(slopeTop / slopeBottom) < 1e-9, fmt("mean %.1e, slope %.1e", mean, slopeTop / slopeBottom));
}