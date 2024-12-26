#include "Log.h"
#include <cstdio> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <filesystem>
#include <chrono>
#include <locale>
#include <vector>
#include "Units.h"


std::string formatWithComma(double value) 
{
    std::ostringstream stream;
    stream.imbue(std::locale("C"));
    stream << std::fixed << std::setprecision(5) << value;
    std::string result = stream.str();

    std::replace(result.begin(), result.end(), '.', ',');
    return result;
}

namespace Log
{
    void avg_sfr(std::vector<Particle*> particles, const double time)
    {
        double totalSFR = 0;
        double n = 0;
        for (size_t i = 0; i < particles.size(); i++)
        {
            if(particles[i]->type == 2) 
            {
                totalSFR += particles[i]->sfr;
                n++;
            }
        }
        Log::printData("total_SFR.csv", time, totalSFR);
        totalSFR /= n;
        Log::printData("avg_SFR.csv", time, totalSFR);
    }
    void total_Mass(std::vector<Particle*> particles, const double time)
    {
        double gasMass = 0;
        double starMass = 0;
        for (size_t i = 0; i < particles.size(); i++)
        {
            if(particles[i]->type == 2) 
            {
                gasMass += particles[i]->mass;
            }
            if(particles[i]->type == 1) starMass += particles[i]->mass;
        }
        Log::printData("gasMass.csv", time, gasMass);
        Log::printData("starMass.csv", time, starMass);
    }
    void avg_U(std::vector<Particle*> particles, const double time)
    {
        double totalU = 0;
        double n = 0;
        for (size_t i = 0; i < particles.size(); i++)
        {
            if(particles[i]->type == 2 && particles[i]->U > 0 && particles[i]->U < 1e11) 
            {
                totalU += particles[i]->U;
                n++;
            }
        }
        Log::printData("total_U.csv", time, totalU);
        totalU /= n;
        Log::printData("avg_U.csv", time, totalU);
    }

    void avg_R_sfr(std::vector<Particle*> particles, int numberOfParticles)
    {
        std::cout << "Saving sfr curve..." << std::endl;
        std::vector<Particle> ps;
        for (int i = 0; i < numberOfParticles; i ++)
        {
            if(particles[i]->type != 2) continue;
            if(!particles[i]) continue;
            Particle p;
            p.position = particles[i]->position;
            p.sfr = particles[i]->sfr;
            ps.push_back(p);
        }
        //get the velocity of the particles
        struct data
        {
            double r;
            double sfr;
            double i;
        };

        std::vector<data> v;
        for (size_t i = 0; i < ps.size(); i++) {
            data d;
            d.r = ps[i].position.length() / Units::KPC;
            d.sfr = ps[i].sfr;
            d.i = i;
            v.push_back(d);
        }

        // Sortiere die Daten nach Abstand
        std::sort(v.begin(), v.end(), [](const data& a, const data& b) {
            return a.r < b.r;
        });

        // Teile die Daten in Bins auf
        double r_abstand = v.back().r / 100;
        std::vector<std::vector<double>> bins(100);

        for (const auto& elem : v) {
            int bin_index = static_cast<int>(elem.r / r_abstand);
            if (bin_index >= 100) bin_index = 99; // Sicherstellen, dass der Index im gültigen Bereich bleibt
            bins[bin_index].push_back(elem.sfr);
        }

        // Berechne den Median der Geschwindigkeiten für jeden Bin
        std::vector<double> medians(100);
        for (size_t i = 0; i < bins.size(); i++) {
            if (bins[i].empty()) {
                medians[i] = 0; // oder einen anderen repräsentativen Wert für leere Bins
                continue;
            }
            size_t mid = bins[i].size() / 2;
            std::nth_element(bins[i].begin(), bins[i].begin() + mid, bins[i].end());
            double median = bins[i][mid];
            if (bins[i].size() % 2 == 0) {
                std::nth_element(bins[i].begin(), bins[i].begin() + mid - 1, bins[i].end());
                median = (median + bins[i][mid - 1]) / 2.0;
            }
            medians[i] = median;
        }
        //print the velocity of the particles
        for (size_t i = 0; i < medians.size(); i++)
        {
            Log::printData("sfr_Curve.csv",(i + 1) * r_abstand, medians[i]);
        }
    }
    void avg_R_U(std::vector<Particle*> particles, int numberOfParticles)
    {
        std::cout << "Saving U curve..." << std::endl;
        std::vector<Particle> ps;
        for (int i = 0; i < numberOfParticles; i ++)
        {
            if(particles[i]->type != 2) continue;
            if(!particles[i]) continue;
            if(particles[i]->U < 0 || particles[i]->U > 1e11) continue;
            Particle p;
            p.position = particles[i]->position;
            p.U = particles[i]->U;
            ps.push_back(p);
        }

        //get the velocity of the particles
        struct data
        {
            double r;
            double U;
            double i;
        };

        std::vector<data> v;
        for (size_t i = 0; i < ps.size(); i++) {
            data d;
            d.r = ps[i].position.length() / Units::KPC;
            d.U = ps[i].U;
            d.i = i;
            v.push_back(d);
        }

        // Sortiere die Daten nach Abstand
        std::sort(v.begin(), v.end(), [](const data& a, const data& b) {
            return a.r < b.r;
        });

        // Teile die Daten in Bins auf
        double r_abstand = v.back().r / 100;
        std::vector<std::vector<double>> bins(100);

        for (const auto& elem : v) {
            int bin_index = static_cast<int>(elem.r / r_abstand);
            if (bin_index >= 100) bin_index = 99; // Sicherstellen, dass der Index im gültigen Bereich bleibt
            bins[bin_index].push_back(elem.U);
        }

        // Berechne den Median der Geschwindigkeiten für jeden Bin
        std::vector<double> medians(100);
        for (size_t i = 0; i < bins.size(); i++) {
            if (bins[i].empty()) {
                medians[i] = 0; // oder einen anderen repräsentativen Wert für leere Bins
                continue;
            }
            size_t mid = bins[i].size() / 2;
            std::nth_element(bins[i].begin(), bins[i].begin() + mid, bins[i].end());
            double median = bins[i][mid];
            if (bins[i].size() % 2 == 0) {
                std::nth_element(bins[i].begin(), bins[i].begin() + mid - 1, bins[i].end());
                median = (median + bins[i][mid - 1]) / 2.0;
            }
            medians[i] = median;
        }
        //print the velocity of the particles
        for (size_t i = 0; i < medians.size(); i++)
        {
            Log::printData("U_Curve.csv",(i + 1) * r_abstand, medians[i]);
        }
    }
    void avg_R_vel(std::vector<Particle*> particles, int numberOfParticles)
    {
        std::cout << "Saving velocity curve..." << std::endl;
        std::vector<Particle> ps;
        for (int i = 0; i < numberOfParticles; i ++)
        {
            if(particles[i]->galaxyPart != 1) continue;
            Particle p;
            p.position = particles[i]->position;
            p.velocity = particles[i]->velocity;
            ps.push_back(p);
        }

        //get the velocity of the particles
        struct data
        {
            double r;
            double v;
            double i;
        };

        std::vector<data> v;
        for (size_t i = 0; i < ps.size(); i++) {
            data d;
            d.r = ps[i].position.length() / Units::KPC;
            d.v = ps[i].velocity.length() / Units::KMS;
            d.i = i;
            v.push_back(d);
        }

        // Sortiere die Daten nach Abstand
        std::sort(v.begin(), v.end(), [](const data& a, const data& b) {
            return a.r < b.r;
        });

        // Teile die Daten in Bins auf
        double r_abstand = v.back().r / 100;
        std::vector<std::vector<double>> bins(100);

        for (const auto& elem : v) {
            int bin_index = static_cast<int>(elem.r / r_abstand);
            if (bin_index >= 100) bin_index = 99; // Sicherstellen, dass der Index im gültigen Bereich bleibt
            bins[bin_index].push_back(elem.v);
        }

        // Berechne den Median der Geschwindigkeiten für jeden Bin
        std::vector<double> medians(100);
        for (size_t i = 0; i < bins.size(); i++) {
            if (bins[i].empty()) {
                medians[i] = 0; // oder einen anderen repräsentativen Wert für leere Bins
                continue;
            }
            size_t mid = bins[i].size() / 2;
            std::nth_element(bins[i].begin(), bins[i].begin() + mid, bins[i].end());
            double median = bins[i][mid];
            if (bins[i].size() % 2 == 0) {
                std::nth_element(bins[i].begin(), bins[i].begin() + mid - 1, bins[i].end());
                median = (median + bins[i][mid - 1]) / 2.0;
            }
            medians[i] = median;
        }
        //print the velocity of the particles
        for (size_t i = 0; i < medians.size(); i++)
        {
            Log::printData("vel_Curve.csv",(i + 1) * r_abstand, medians[i]);
        }
    }


    std::string outputDir;
    void setOutputDir(const std::string& dir) 
    {
        outputDir = dir;
        if (!std::filesystem::exists(outputDir)) 
        {
            std::filesystem::create_directories(outputDir);
        }
    }

//save data in csv file

    void printData(const std::string& filename, const double x, const double y) 
    { 
        std::ofstream file(outputDir + "/"+ filename, std::ios::app);

        // Überprüfen, ob die Datei erfolgreich geöffnet wurde
        if (!file.is_open()) {
            throw std::ios_base::failure("Fehler beim Öffnen der Datei.");
        }

        // Formatieren der Zahlen mit Punkt als Dezimaltrennzeichen
        std::ostringstream formattedX, formattedY;
        formattedX << std::fixed << std::setprecision(10) << x;
        formattedY << std::fixed << std::setprecision(10) << y;
        // Ersetzen des Punktes durch ein Komma
        std::string xStr = formattedX.str();
        std::string yStr = formattedY.str();
        std::replace(xStr.begin(), xStr.end(), '.', ',');
        std::replace(yStr.begin(), yStr.end(), '.', ',');

        // Schreiben der formatierten Daten in die Datei
        file << xStr << ";" << yStr << "\n";

        // Datei schließen
        file.close();
    }
    
//track process time
    bool hasStarted = false;
    std::ofstream LogsDir;
    std::ofstream proccessFile;
    std::chrono::steady_clock::time_point startTimestamp;
    std::string currentProcessName;


    void startProcess(const std::string& processName) 
    {
        if(hasStarted == true)
        {
            endProcess();
        }
        else
        {
            std::string filename =  outputDir + "/processLog.csv";

            if (std::ifstream(filename)) 
            {
                std::cout << "Logfile exists. Deleting existing logfile: " << filename << std::endl;
                if (std::remove(filename.c_str()) != 0) 
                {
                    std::cerr << "Error deleting logfile." << std::endl;
                }
            }

            proccessFile.open(filename, std::ios::out | std::ios::app);

            if (proccessFile.is_open()) 
            {
                //std::cout << "Logfile opened successfully." << std::endl;
            } 
            else 
            {
                std::cerr << "Error opening logfile." << std::endl;
        }
            hasStarted = true;
        }
        currentProcessName = processName;
        startTimestamp = std::chrono::steady_clock::now();
    }

    void endProcess() 
    {
        if (!currentProcessName.empty()) {
            auto now = std::chrono::steady_clock::now();
            double elapsedTime = std::chrono::duration<double>(now - startTimestamp).count();
            proccessFile << currentProcessName << ";" << formatWithComma(elapsedTime) << "\n";
            
            proccessFile.flush();

            currentProcessName.clear();
        }
    }
}