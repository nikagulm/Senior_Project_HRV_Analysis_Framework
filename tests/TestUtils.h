#ifndef TESTUTILS_H
#define TESTUTILS_H
#include <cstdio>
#include <string>
inline int passed = 0;
inline int failed = 0;

inline void check(const std::string& name, bool ok, const std::string& detail) {
    printf("%-4s  %-58s %s\n", ok ? "PASS" : "FAIL", name.c_str(), detail.c_str());
    if (ok) {
        passed++;
    } else {
        failed++;
    }
}

inline std::string fmt(const char* format, double a, double b = 0, double c = 0) {
    char text[200];
    snprintf(text, sizeof(text), format, a, b, c);
    return text;
}

void testTransforms();
void testGoertzel();
void testBandPowers();
void testTimeDomain();
void testNNRules();
void testSpline();
void testEvaluator();
void testDetectors();
void testRecords();
void testWrongInput();

#endif