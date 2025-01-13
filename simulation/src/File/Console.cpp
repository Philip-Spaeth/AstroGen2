#include "Console.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/ioctl.h>

using namespace std;
namespace fs = std::filesystem;

void Console::printProgress(double currentStep, double steps, std::string text) 
{
    // Verhindert das Aktualisieren, wenn der aktuelle Schritt größer als die Gesamt-Schritte ist
    if(currentStep > steps) return;

    // Initialisierung des Timers beim ersten Aufruf
    if (!timerStarted) {
        std::cout << std::endl;
        startTime = std::chrono::high_resolution_clock::now();
        timerStarted = true;
    }

    // Berechnung des Fortschritts
    double progress = currentStep / steps;

    // Terminalgröße ermitteln
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    int terminalWidth = w.ws_col;

    // Definieren der statischen Teile der Ausgabe
    std::string prefix = "Progress:";
    std::string suffix = text;
    std::string fixedParts = " []  %  Estimated remaining time:  "; // Platzhalter für Balken, Prozent und Zeit

    // Maximale Anzahl der Zeichen für die Progress-Bar berechnen
    // Reserviert Platz für Prefix, Suffix und feste Teile der Ausgabe
    int reservedSpace = prefix.length() + suffix.length() + fixedParts.length() + 10; // 10 für Prozentzahl und Zeitangabe
    int barWidth = terminalWidth - reservedSpace;

    if (barWidth < 10) barWidth = 10; // Mindestbreite der Progress-Bar

    // Berechnung der Position der Fortschrittsanzeige im Balken
    int pos = static_cast<int>(barWidth * progress);

    // Aufbau der Progress-Bar
    std::ostringstream barStream;
    barStream << prefix << " [";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) {
            barStream << "=";
        } else if (i == pos) {
            barStream << ">";
        } else {
            barStream << " ";
        }
    }
    barStream << "] ";

    // Berechnung des Prozentsatzes
    int percent = static_cast<int>(progress * 100.0);

    // Berechnung der verbleibenden Zeit
    double remainingTime = 0.0;

    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = currentTime - startTime;

    if (currentStep > 0) { // Vermeidung von Division durch Null
        double estimatedTotalTime = (elapsed.count() / currentStep) * steps;
        double estimatedRemainingTime = estimatedTotalTime - elapsed.count();
        remainingTime = estimatedRemainingTime;
    }

    // Formatierung der verbleibenden Zeit
    std::string timeUnit;
    if (remainingTime < 60) {
        timeUnit = "s";
    } else if (remainingTime < 3600) {
        remainingTime /= 60;
        timeUnit = "min";
    } else {
        remainingTime /= 3600;
        timeUnit = "h";
    }

    // Formatierung der verbleibenden Zeit mit einer Nachkommastelle
    std::ostringstream timeStream;
    timeStream << std::fixed << std::setprecision(1) << remainingTime;
    std::string timeleft = timeStream.str();

    // Aufbau der gesamten Ausgabe
    barStream << percent << "%  Estimated remaining time: " << timeleft << timeUnit << "  " << suffix;

    // Sicherstellen, dass die gesamte Ausgabe die Terminalbreite nicht überschreitet
    std::string progressBarStr = barStream.str();
    if (progressBarStr.length() > static_cast<size_t>(terminalWidth)) {
        // Kürzen der Suffix- oder anderer Teile, falls nötig
        int excess = progressBarStr.length() - terminalWidth;
        if (excess > 0 && suffix.length() > static_cast<size_t>(excess)) {
            barStream.str("");
            barStream << prefix << " [";
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) {
                    barStream << "=";
                } else if (i == pos) {
                    barStream << ">";
                } else {
                    barStream << " ";
                }
            }
            barStream << "] ";

            barStream << percent << "%  Estimated remaining time: " << timeleft << timeUnit << "  " 
                      << suffix.substr(0, suffix.length() - excess) << "...";
            progressBarStr = barStream.str();
        }
    }

    // Ausgabe der Progress-Bar mit Carriage Return zur Aktualisierung der gleichen Zeile
    std::cout << "\r" << progressBarStr << std::flush;

    // Beenden der Progress-Bar nach dem letzten Schritt
    if (currentStep >= steps) {
        std::cout << std::endl;
    }
}

void Console::printSystemInfo()
{
    std::cout << "\nComputational parameters:" << std::endl;
    // Number of CPU cores
    int numCores = 0;
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    numCores = sysinfo.dwNumberOfProcessors;
    #else
    numCores = sysconf(_SC_NPROCESSORS_ONLN);
    #endif

    std::cout << "  Number of CPU Cores: " << numCores << std::endl;

    // CPU Clock Speed
    double cpuFrequency = 0.0;
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.find("cpu MHz") != std::string::npos) {
            cpuFrequency = std::stod(line.substr(line.find(":") + 1));
            break;
        }
    }

    //display the frequency in GHz iwth 2 decimal places
    std::cout << "  CPU Clock Speed: " << std::fixed << std::setprecision(2) << cpuFrequency / 1000 << " GHz" << std::endl;

    // Available RAM
    size_t totalRam = 0;
    std::ifstream meminfo("/proc/meminfo");
    while (std::getline(meminfo, line)) {
        if (line.find("MemTotal") != std::string::npos) {
            totalRam = std::stoul(line.substr(line.find(":") + 1)) / 1024; // in MB
            break;
        }
    }

    std::cout << "  Available RAM: " << std::fixed << std::setprecision(1) << (double)totalRam / 1000.0 << " GB" << std::endl;

    // Storage Write Speed
    const size_t fileSize = 100 * 1024 * 1024; // 100 MB
    char *buffer = new char[fileSize];
    std::fill(buffer, buffer + fileSize, 'A');

    auto start = std::chrono::high_resolution_clock::now();

    std::ofstream file("test_file.bin", std::ios::binary);
    file.write(buffer, fileSize);
    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    double writeSpeed = fileSize / duration.count() / (1024 * 1024); // MB/s

    std::cout << "  Storage Write Speed: " << writeSpeed << " MB/s \n" << std::endl;

    delete[] buffer;
}