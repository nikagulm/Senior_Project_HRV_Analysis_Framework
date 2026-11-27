#include "TestUtils.h"

int main() {
    testTransforms();
    testGoertzel();
    testBandPowers();
    testTimeDomain();
    testNNRules();
    testSpline();
    testEvaluator();
    testDetectors();
    testRecords();
    testWrongInput();
    printf("\n%d tests passed, %d failed\n", passed, failed);
    return failed == 0 ? 0 : 1;
}