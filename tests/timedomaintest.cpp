#include <iostream>
#include <vector>
#include "../TimeDomainHRV.h"
using namespace std;

int main() {
    vector<double> nn = {0.80, 0.85, 0.78, 0.90, 0.82};
    TimeDomainHRV td;
    TimeDomainResult r = td.calculate(nn);
    cout << "meanNN= " << r.meanNN << " ms (expected 830)" << endl;
    cout << "sdnn= " << r.sdnn   << " ms (expected 46.9)" << endl;
    cout << "RMSSD= " << r.rmssd  << " ms (expected 84.0)" << endl;
    cout << "pNN50= " << r.pnn50  << " % (expected 75)" << endl;
}