#include "annotationParser.h"
#include <fstream>
#include <iostream>

using namespace std;


// Checks whether the annotation type represents an actual heartbeat
bool isBeatAnnotation(int type) {
    return (type >= 1 && type <= 13) ||
           type == 25 ||
           type == 30 ||
           type == 34 ||
           type == 35 ||
           type == 38 ||
           type == 41;
}


/**
 * Reads reference heartbeat annotations from a PhysioNet .atr file.
 *
 * @param fileName Path to the annotation file.
 * @return Sample positions of the reference heartbeats.
 */
vector<int> AnnotationParser::parse(const string& fileName) {
    ifstream file(fileName, ios::binary);

    if (!file.is_open()) {
        cout << "Can't open the annotation file" << endl;
        return {};
    }

    vector<int> annotations;
    unsigned char bytes[2];

    int currentSample = 0;

    while (file.read(reinterpret_cast<char*>(bytes), 2)) {

        // Lower 10 bits contain the time difference
        int timeDifference = bytes[0] | ((bytes[1] & 0x03) << 8);

        // Upper 6 bits contain the annotation type
        int annotationType = bytes[1] >> 2;


        // SKIP: a larger time offset is stored in the next 4 bytes
        if (annotationType == 59) {
            unsigned char skipBytes[4];

            if (!file.read(reinterpret_cast<char*>(skipBytes), 4)) {
                break;
            }

            // Combine the four bytes into the larger time offset
            int skip = (skipBytes[0] << 16) | (skipBytes[1] << 24)
                | skipBytes[2] | (skipBytes[3] << 8);

            currentSample += skip;
            continue;
        }


        // End of annotation file
        if (annotationType == 0 && timeDifference == 0) {
            break;
        }


        // NUM, SUB and CHN contain metadata rather than beat positions
        if (annotationType == 60 ||
            annotationType == 61 ||
            annotationType == 62) {
            continue;
        }


        // AUX is followed by additional information
        if (annotationType == 63) {
            int auxLength = timeDifference;

            // AUX data is padded to an even number of bytes
            if (auxLength % 2 != 0) {
                auxLength++;
            }

            file.ignore(auxLength);
            continue;
        }


        // Convert the relative time to an absolute sample position
        currentSample += timeDifference;


        // Keep only actual heartbeat annotations
        if (isBeatAnnotation(annotationType)) {
            annotations.push_back(currentSample);
        }
    }

    return annotations;
}