#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <memory>
#include "Particle.h"



namespace Log
{
    void setOutputDir(const std::string& outputDir);
    extern std::string outputDir;

//logs are saved in the output dir /logs
    // track process time
    void startProcess(const std::string& processName);
    void endProcess();
    extern bool hasStarted;
    extern std::ofstream LogsDir;
    extern std::ofstream proccessFile;
    extern std::chrono::steady_clock::time_point startTimestamp;
    extern std::string currentProcessName;

    //save data to in csv file
    void printData(const std::string& file, const double x, const double y);

    //save velocity for each radius, rotation curve
    void avg_R_vel(std::vector<Particle*> particles, int numberOfParticles);

    //save SFR for each radius
    void avg_R_sfr(std::vector<Particle*> particles, int numberOfParticles);
    //save avrage SFR for the whole system
    void sfr(std::vector<Particle*> particles, const double time, double sfr);
    
    //save total mass for the whole system
    void total_Mass(std::vector<Particle*> particles, const double time);
    
    //save U for each radius
    void avg_R_U(std::vector<Particle*> particles, int numberOfParticles);
    //save avrage U for the whole system
    void avg_U(std::vector<Particle*> particles, const double time);
}
#endif