#ifndef DATAPARSER_H
#define DATAPARSER_H
#include <string>
#include <vector>

class DataParser {
public:
    std::vector<int> parseFormat212(const std::string& fileName);
    std::vector<int> getChannel1() const;
    std::vector<int> getChannel2() const;

private:
    std::vector<int> channel1;
    std::vector<int> channel2;

};

#endif