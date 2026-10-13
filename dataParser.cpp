#include "dataParser.h"
#include <fstream>
#include <iostream>
using namespace std;


/* reads a PhysioNet .dat file in Format 212
   decodes two signed 12-bit ECG samples from every 3 bytes and returns them as integer values
   - @param fileName Path to the .dat file.
   - @return Decoded ECG samples.
*/
vector<int> DataParser::parseFormat212(const string& fileName) {
    ifstream file(fileName, ios::binary);
    if (!file.is_open()) {
        cout << "Can't open the data file" << endl;
        return {};
    }

    vector<int> samples;

    channel1.clear();
    channel2.clear();
    unsigned char bytes[3];
    while (file.read(reinterpret_cast<char*>(bytes), 3)) {
        // decode two 12-bit samples from these 3 bytes here.
        int sample1 = bytes[0] | ((bytes[1] & 0x0F) << 8);
        int sample2 = bytes[2] | ((bytes[1] & 0xF0) << 4);
        if (sample1 >= 2048) {
            sample1 -= 4096;
            }

        if (sample2 >= 2048) {
            sample2 -= 4096;    
            }
        
        samples.push_back(sample1);
        samples.push_back(sample2);

        channel1.push_back(sample1);
        channel2.push_back(sample2);
    }
    return samples;
}

vector<int> DataParser::getChannel1() const {
    return channel1;
}

vector<int> DataParser::getChannel2() const {
    return channel2;
}