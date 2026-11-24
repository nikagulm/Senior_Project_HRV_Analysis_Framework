#ifndef REPORTWRITER_H
#define REPORTWRITER_H
#include <string>
#include <vector>
#include "detectorEvaluator.h"
#include "TimeDomainHRV.h"
#include "FrequencyDomainHRV.h"

struct ReportRow {
    std::string record;
    std::string detector;
    bool isReference; 
    EvaluationResult eval;
    TimeDomainResult time;
    FrequencyDomainResult freq;
    double timeMs; 
};

struct PlotSeries {
    std::string record;
    std::string detector;
    std::vector<double> peakTimes; 
    std::vector<double> nnTimes;  
    std::vector<double> nn;        
};

class ReportWriter {
public:
    void addReference(const std::string& record, const TimeDomainResult& time, const FrequencyDomainResult& freq);
    void addDetector(const std::string& record, const std::string& detector, const EvaluationResult& eval, const TimeDomainResult& time, const FrequencyDomainResult& freq, double timeMs);
    void printRecord(const std::string& record);
    void printSummary();
    bool writeCSV(const std::string& fileName);
    void writeECG(const std::string& folder, const std::string& record, const std::vector<double>& ecg, int samplingRate);
    void addPlotSeries(const PlotSeries& series);
    bool writePlotData(const std::string& folder);

private:
    std::vector<ReportRow> rows;
    std::vector<PlotSeries> series;
    std::vector<std::string> ecgRecords;
    std::vector<int> ecgRates;
    const ReportRow* findReference(const std::string& record);
    double percentDifference(double value, double reference);
};

#endif