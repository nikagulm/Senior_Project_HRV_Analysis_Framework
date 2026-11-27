//reading the physionet records and handling wrong or missing input
#include "TestUtils.h"
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <filesystem>
#include "headerParser.h"
#include "dataParser.h"
#include "annotationParser.h"
using namespace std;

//physionet checksum ->> sum of all samples of a signal kept as a 16-bit number
int checksum(const vector<int>& samples) {
    long long sum = 0;
    for (int s : samples) {
        sum += s;
    }
    return static_cast<int16_t>(sum & 0xFFFF);
}

void testRecords() {
    printf("\nReading the PhysioNet records\n");
    string folder = "records_physionet/";
    if (!filesystem::exists(folder + "100.hea")) {
        printf("records_physionet/ not found, record tests skipped (run from the project folder)\n");
        return;
    }
    HeaderParser h;
    h.parse(folder + "100.hea");
    check("header of record 100 (360 Hz, 650000 samples, 2 signals)", h.getSamplingRate() == 360 && h.getNumberOfSamples() == 650000 && h.getSignals().size() == 2, fmt("%.0f Hz, %.0f samples, %.0f signals", h.getSamplingRate(), h.getNumberOfSamples(), h.getSignals().size()));
    DataParser d;
    d.parseFormat212(folder + "100.dat");
    const SignalInfo& s0 = h.getSignals()[0];
    check("record 100: first sample equals the initial value in the header", !d.getChannel1().empty() && d.getChannel1()[0] == s0.initialValue, fmt("decoded %.0f, header %.0f", d.getChannel1().empty() ? 0 : d.getChannel1()[0], s0.initialValue));
    AnnotationParser a;
    vector<int> beats = a.parse(folder + "100.atr");
    check("record 100: number of annotated beats (expected 2273)", beats.size() == 2273, fmt("%.0f beats", beats.size()));
    vector<string> records = {
        "100", "101", "102", "103", "104", "105", "106", "107", "108", "109",
        "111", "112", "113", "114", "115", "116", "117", "118", "119",
        "121", "122", "123", "124", "200", "201", "202", "203", "205", "207", "208", "209", "210",
        "212", "213", "214", "215", "217", "219", "220", "221", "222", "223",
        "228", "230", "231", "232", "233", "234", "118e24", "118e18", "118e12", "118e06", "118e00", "118e_6",
        "119e24", "119e18", "119e12", "119e06", "119e00", "119e_6"
    };
    int ok = 0, format16 = 0;
    string bad;
    for (const string& r : records) {
        HeaderParser hr;
        hr.parse(folder + r + ".hea");
        if (hr.getSignals().size() < 2) {
            bad += r + " ";
            continue;
        }
        DataParser dr;
        if (hr.getSignals()[0].format == 16) {
            dr.parseFormat16(folder + r + ".dat");
            format16++;
        } else {
            dr.parseFormat212(folder + r + ".dat");
        }
        bool good = dr.getChannel1().size() == hr.getNumberOfSamples() && checksum(dr.getChannel1()) == hr.getSignals()[0].checksum && checksum(dr.getChannel2()) == hr.getSignals()[1].checksum;
        if (good) {
            ok++;
        } else {
            bad += r + " ";
        }
    }
    check("sample count and checksums of both channels, all records", ok == records.size(), fmt("%.0f of %.0f records (%.0f in format 16)", ok, records.size(), format16) + (bad.empty() ? "" : ", failed: " + bad));
}

//wrong or missing input
void testWrongInput() {
    printf("\nWrong or missing input\n");
    string tmp = filesystem::temp_directory_path().string() + "/hrv_test_";
    //missing files
    HeaderParser h1;
    h1.parse(tmp + "missing.hea");
    check("missing header file: no signals, no crash", h1.getSignals().empty(), "");
    DataParser d1;
    check("missing signal file: no samples, no crash", d1.parseFormat212(tmp + "missing.dat").empty(), "");
    AnnotationParser a1;
    check("missing annotation file: no beats, no crash", a1.parse(tmp + "missing.atr").empty(), "");

    //empty header file
    ofstream(tmp + "empty.hea").close();
    HeaderParser h2;
    h2.parse(tmp + "empty.hea");
    check("empty header file: no signals, no crash", h2.getSignals().empty(), "");

    //a web page saved instead of the record (happened during the download)
    ofstream(tmp + "page.hea") << "<html>\n<head><title>404 Not Found</title></head>\n<body>Not Found</body>\n</html>\n";
    bool crashed = false;
    HeaderParser h3;
    try {
        h3.parse(tmp + "page.hea");
    } catch (...) {
        crashed = true;
    }
    check("web page instead of a header: no signals, no crash", !crashed && h3.getSignals().empty(), crashed ? "exception thrown" : "");
    //incomplete signal file
    string folder = "records_physionet/";
    if (filesystem::exists(folder + "100.dat")) {
        ifstream in(folder + "100.dat", ios::binary);
        vector<char> part(300000);
        in.read(part.data(), part.size());
        ofstream(tmp + "cut.dat", ios::binary).write(part.data(), part.size());
        DataParser d2;
        d2.parseFormat212(tmp + "cut.dat");
        check("cut-off signal file is recognised by the sample count", d2.getChannel1().size() != 650000, fmt("%.0f samples instead of 650000", d2.getChannel1().size()));
    }
}