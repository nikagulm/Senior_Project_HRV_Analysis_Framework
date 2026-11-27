//matching detections with reference beats
#include "TestUtils.h"
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include "detectorEvaluator.h"
using namespace std;

EvaluationResult slowEvaluate(const vector<int>& det, const vector<int>& ann, int tol, int start) {
    EvaluationResult r{};
    vector<bool> used(ann.size(), false);
    for (int d : det) {
        int best = -1;
        for (int i = 0; i < ann.size(); i++) {
            if (!used[i] && abs(d - ann[i]) <= tol && (best == -1 || abs(d - ann[i]) < abs(d - ann[best]))) {
                best = i;
            }
        }
        if (best != -1) {
            used[best] = true;
            if (ann[best] >= start) {
                r.truePositives++;
            }
        } else if (d >= start) {
            r.falsePositives++;
        }
    }
    for (int i = 0; i < ann.size(); i++) {
        if (!used[i] && ann[i] >= start) {
            r.falseNegatives++;
        }
    }
    return r;
}

void testEvaluator() {
    printf("\nDetector evaluation\n");
    DetectorEvaluator ev;
    //1010 matches 1000, 2100 is too far from 2000, 3000 matches, 3500 has no annotation
    EvaluationResult r = ev.evaluate({1010, 2100, 3000, 3500}, {1000, 2000, 3000}, 54, 0);
    check("hand example (expected TP 2, FP 2, FN 1)", r.truePositives == 2 && r.falsePositives == 2 && r.falseNegatives == 1, fmt("TP %.0f, FP %.0f, FN %.0f", r.truePositives, r.falsePositives, r.falseNegatives));
    check("  sensitivity 2/3 and PPV 2/4", abs(r.sensitivity - 2.0 / 3) < 1e-12 && abs(r.ppv - 0.5) < 1e-12, fmt("Se %.4f, PPV %.4f", r.sensitivity, r.ppv));

    //beats before the start sample are matched but not counted
    EvaluationResult s = ev.evaluate({500, 1000, 2000}, {500, 1000, 2000}, 54, 900);
    check("beats before the start sample are not counted", s.truePositives == 2 && s.falsePositives == 0, fmt("TP %.0f, FP %.0f", s.truePositives, s.falsePositives));

    //linear version must agree with the straightforward one on random cases
    mt19937 rng(3);
    int same = 0, cases = 20000;
    for (int c = 0; c < cases; c++) {
        vector<int> ann, det;
        int t = 0;
        int beats = 5 + rng() % 40;
        for (int i = 0; i < beats; i++) {
            t += 72 + rng() % 400; //at least 200 ms between beats
            ann.push_back(t);
            if (rng() % 10) {
                det.push_back(t + static_cast<int>(rng() % 121) - 60); //sometimes outside the tolerance
            }
            if (rng() % 8 == 0) {
                det.push_back(t + 100 + rng() % 200); //false beat
            }
        }
        sort(det.begin(), det.end());
        int start = rng() % (t / 2 + 1);
        EvaluationResult a = ev.evaluate(det, ann, 54, start);
        EvaluationResult b = slowEvaluate(det, ann, 54, start);
        if (a.truePositives == b.truePositives && a.falsePositives == b.falsePositives && a.falseNegatives == b.falseNegatives) {
            same++;
        }
    }
    check("linear evaluator equals the quadratic one on random cases", same == cases, fmt("%.0f of %.0f identical", same, cases));
}