#ifndef TIMER_H
#define TIMER_H
#include <chrono>

class Timer {
public:
    void start() {
        begin = std::chrono::steady_clock::now();
    }
    double stop() {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        return std::chrono::duration<double, std::milli>(end - begin).count();
    }

private:
    std::chrono::steady_clock::time_point begin;
};

#endif