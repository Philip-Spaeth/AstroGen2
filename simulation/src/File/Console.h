#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <chrono>

class Console
{
public:
    Console(){};
    ~Console(){}
///console output
    //progress bar
    void printProgress(double currentStep, double steps, std::string text);
    //system info
    static void printSystemInfo();

private:
    std::chrono::_V2::system_clock::time_point startTime;
    bool timerStarted = false;
};