#include "detectorEvaluator.h"
#include <cmath>
using namespace std;

EvaluationResult DetectorEvaluator::evaluate(
    const vector<int>& detectedPeaks,
    const vector<int>& annotations,
    int toleranceSamples,
    int startSample

) {
    EvaluationResult result{};
    vector<bool> matchedAnnotations(annotations.size(), false);

    //check every detected peak
    for (int detected : detectedPeaks) {
        int bestMatch = -1;
        int smallestDifference = toleranceSamples + 1;

        //find the closest unmatched reference annotation
        for (int i = 0; i < annotations.size(); i++) {
            if (matchedAnnotations[i]) {
                continue;
            }
            int difference = abs(detected - annotations[i]);

            //keep the closest annotation within the allowed tolerance
            if (difference <= toleranceSamples &&
                difference < smallestDifference) {
                smallestDifference = difference;
                bestMatch = i;
            }
        }

        //if a match was found, count a true positive
        if (bestMatch != -1) {
            matchedAnnotations[bestMatch] = true;
            if (annotations[bestMatch] >= startSample){
                result.truePositives++;
            }
        } else if (detected >= startSample) {
            result.falsePositives++;
        }
    }

    //count reference beats that were never matched after first 5 mins
     for (int i = 0; i < annotations.size(); i++) {
        if (!matchedAnnotations[i] && annotations[i] >= startSample) {
            result.falseNegatives++;
         }
     }

    if (result.truePositives + result.falseNegatives > 0) {
        result.sensitivity = static_cast<double>(result.truePositives) / (result.truePositives + result.falseNegatives);
    }

    if (result.truePositives + result.falsePositives > 0) {
        result.ppv = static_cast<double>(result.truePositives) / (result.truePositives + result.falsePositives);
    }

    return result;
}
