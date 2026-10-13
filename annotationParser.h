#ifndef ANNOTATIONPARSER_H
#define ANNOTATIONPARSER_H

#include <string>
#include <vector>

class AnnotationParser {
public:
    std::vector<int> parse(const std::string& fileName);
};

#endif