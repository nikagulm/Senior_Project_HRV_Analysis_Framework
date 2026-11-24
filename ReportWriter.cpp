#include "ReportWriter.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <filesystem>
using namespace std;

void ReportWriter::addReference(const string& record, const TimeDomainResult& time, const FrequencyDomainResult& freq) {
    ReportRow row;
    row.record = record;
    row.detector = "Reference";
    row.isReference = true;
    row.eval = {0, 0, 0, 0.0, 0.0};
    row.time = time;
    row.freq = freq;
    row.timeMs = 0.0;
    rows.push_back(row);
}

void ReportWriter::addDetector(const string& record, const string& detector, const EvaluationResult& eval, const TimeDomainResult& time, const FrequencyDomainResult& freq, double timeMs) {
    ReportRow row;
    row.record = record;
    row.detector = detector;
    row.isReference = false;
    row.eval = eval;
    row.time = time;
    row.freq = freq;
    row.timeMs = timeMs;
    rows.push_back(row);
}

const ReportRow* ReportWriter::findReference(const string& record) {
    for (const ReportRow& row : rows) {
        if (row.record == record && row.isReference) {
            return &row;
        }
    }
    return nullptr;
}

//how far a value is from the reference in %
double ReportWriter::percentDifference(double value, double reference) {
    if (reference == 0.0) {
        return 0.0;
    }
    return 100.0 * (value - reference) / reference;
}

void ReportWriter::printRecord(const string& record) {
    const ReportRow* ref = findReference(record);

    printf("\nRecord %s\n", record.c_str());
    printf("%-20s %7s %7s %4s %4s | %7s %6s %6s %6s | %7s %7s %6s | %7s %7s | %8s\n",
           "Detector", "Se%", "PPV%", "FP", "FN", "meanNN", "SDNN", "RMSSD", "pNN50",
           "LF", "HF", "LF/HF", "dRMSSD%", "dLF%", "Time ms");
    printf("%s\n", string(137, '-').c_str());

    for (const ReportRow& row : rows) {
        if (row.record != record) {
            continue;
        }

        //detection scores only exist for detectors
        if (row.isReference) {
            printf("%-20s %7s %7s %4s %4s", row.detector.c_str(), "-", "-", "-", "-");
        } else {
            printf("%-20s %7.2f %7.2f %4d %4d", row.detector.c_str(), row.eval.sensitivity * 100, row.eval.ppv * 100, row.eval.falsePositives, row.eval.falseNegatives);
        }

        printf(" | %7.1f %6.1f %6.1f %6.2f | %7.1f %7.1f %6.3f", row.time.meanNN, row.time.sdnn, row.time.rmssd, row.time.pnn50, row.freq.lf, row.freq.hf, row.freq.lfhf);

        //difference from the reference HRV
        if (row.isReference || ref == nullptr) {
            printf(" | %7s %7s", "-", "-");
        } else {
            printf(" | %+7.1f %+7.1f", percentDifference(row.time.rmssd, ref->time.rmssd), percentDifference(row.freq.lf, ref->freq.lf));
        }

        //detection time
        if (row.isReference) {
            printf(" | %8s\n", "-");
        } else {
            printf(" | %8.1f\n", row.timeMs);
        }
    }
}

void ReportWriter::printSummary() {
    //detector names in the order they first appear
    vector<string> detectors;
    for (const ReportRow& row : rows) {
        if (row.isReference) {
            continue;
        }
        bool seen = false;
        for (const string& d : detectors) {
            if (d == row.detector) {
                seen = true;
            }
        }
        if (!seen) {
            detectors.push_back(row.detector);
        }
    }

    printf("\nSummary over all records (average; HRV errors are absolute %% differences from the reference)\n");
    printf("%-20s %7s %7s %6s %6s | %9s %9s %9s | %8s\n", "Detector", "Se%", "PPV%", "FP", "FN", "|dRMSSD|%", "|dLF|%", "|dHF|%", "Time ms");
    printf("%s\n", string(97, '-').c_str());

    for (const string& d : detectors) {
        double se = 0, ppv = 0, rmssdError = 0, lfError = 0, hfError = 0, time = 0;
        int fp = 0, fn = 0, count = 0;
        for (const ReportRow& row : rows) {
            if (row.detector != d) {
                continue;
            }
            const ReportRow* ref = findReference(row.record);
            se += row.eval.sensitivity * 100;
            ppv += row.eval.ppv * 100;
            fp += row.eval.falsePositives;
            fn += row.eval.falseNegatives;
            if (ref != nullptr) {
                rmssdError += fabs(percentDifference(row.time.rmssd, ref->time.rmssd));
                lfError += fabs(percentDifference(row.freq.lf, ref->freq.lf));
                hfError += fabs(percentDifference(row.freq.hf, ref->freq.hf));
            }
            time += row.timeMs;
            count++;
        }

        printf("%-20s %7.2f %7.2f %6d %6d | %9.1f %9.1f %9.1f | %8.1f\n", d.c_str(), se / count, ppv / count, fp, fn, rmssdError / count, lfError / count, hfError / count, time / count);
    }
}

bool ReportWriter::writeCSV(const string& fileName) {
    //make the folder if it does not exist yet
    filesystem::path path(fileName);
    if (path.has_parent_path()) {
        filesystem::create_directories(path.parent_path());
    }

    ofstream file(fileName);
    if (!file.is_open()) {
        printf("Can't write %s\n", fileName.c_str());
        return false;
    }

    file << "record,detector,TP,FP,FN,Se,PPV,meanNN,SDNN,RMSSD,pNN50,VLF,LF,HF,LF_HF,segments,time_ms\n";
    for (const ReportRow& row : rows) {
        file << row.record << "," << row.detector << ",";
        if (row.isReference) {
            file << ",,,,,";
        } else {
            file << row.eval.truePositives << "," << row.eval.falsePositives << "," << row.eval.falseNegatives << "," << row.eval.sensitivity * 100 << "," << row.eval.ppv * 100 << ",";
        }
        file << row.time.meanNN << "," << row.time.sdnn << "," << row.time.rmssd << "," << row.time.pnn50 << "," << row.freq.vlf << "," << row.freq.lf << "," << row.freq.hf << "," << row.freq.lfhf << "," << row.freq.segments << ",";
        if (!row.isReference) {
            file << row.timeMs;
        }
        file << "\n";
    }
    printf("\nResults saved to %s\n", fileName.c_str());
    return true;
}

//whole ECG as raw float32, small and fast to load in the plot script
void ReportWriter::writeECG(const string& folder, const string& record, const vector<double>& ecg, int samplingRate) {
    filesystem::create_directories(folder);
    ofstream file(folder + "/ecg_" + record + ".bin", ios::binary);
    for (double x : ecg) {
        float value = static_cast<float>(x);
        file.write(reinterpret_cast<const char*>(&value), sizeof(float));
    }
    ecgRecords.push_back(record);
    ecgRates.push_back(samplingRate);
}

void ReportWriter::addPlotSeries(const PlotSeries& s) {
    series.push_back(s);
}

//one CSV per kind of data each with a record and detector column
bool ReportWriter::writePlotData(const string& folder) {
    filesystem::create_directories(folder);

    ofstream infoFile(folder + "/records.csv");
    ofstream peakFile(folder + "/peaks.csv");
    ofstream nnFile(folder + "/nn.csv");
    if (!infoFile.is_open() || !peakFile.is_open() || !nnFile.is_open()) {
        printf("Can't write the plot data to %s\n", folder.c_str());
        return false;
    }

    infoFile << "record,fs\n";
    for (int r = 0; r < ecgRecords.size(); r++) {
        infoFile << ecgRecords[r] << "," << ecgRates[r] << "\n";
    }

    peakFile << "record,detector,time\n";
    nnFile << "record,detector,time,nn_ms\n";
    for (const PlotSeries& s : series) {
        for (double t : s.peakTimes) {
            peakFile << s.record << "," << s.detector << "," << t << "\n";
        }
        for (int i = 0; i < s.nn.size(); i++) {
            nnFile << s.record << "," << s.detector << "," << s.nnTimes[i] << "," << s.nn[i] * 1000.0 << "\n";
        }
    }
    printf("Plot data saved to %s\n", folder.c_str());
    return true;
}