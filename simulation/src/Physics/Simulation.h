#pragma once
#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>
#include "Constants.h"
#include "TimeIntegration.h"
#include "DataManager.h"
#include "Console.h"
#include "random.h"
#include "Tree/Node.h"
#include "Tree/Tree.h"
#include <thread>
#include <mutex>
#include <atomic>
#include "Cooling.h"
#include "SFR.h"
#include <chrono>

class Tree;
class DataManager;
class Halo;

class Simulation
{
public:
    Simulation();
    ~Simulation();
    bool init();
    void run();

    //simulation parameters, has to the same as in the input dataset(ICs)
    int numberOfParticles;

    int numParticlesOutput; // Number of particles to output and save in the data files

    //adaptive time integration
    double eta;      // Accuracy parameter for adaptive time step
    double maxTimeStep; // Maximum allowed time step
    double minTimeStep; // Minimum allowed time step

    double globalTime; // global time of the simulation in s
    double endTime; //end time of the simulation in s

    std::chrono::high_resolution_clock::time_point startTimeSimulation; //start time of the simulation

    //save data at each maxTimeStep
    double fixedTimeSteps; //number of fixed maxtime steps
    double fixedStep; //time step in s

    //gravitational softening, adapt it to the size of the system, osftening beginns at 2.8 * e0
    double e0; //softening factor
    
    //SPH parameters
    double massInH; //in kg

    bool starFormation = false;
    bool coolingEnabled = false;

    //Visual density, for all particles, just for visualization, beacuse the real density is only calculated for Gas particles, has no physical meaning
    double visualDensityRadius; //in m

    //Constant hubble expansion
    double H0; //Hubble constant in km/s/Mpc

    //octree with all particles
    double theta;
    
    //particles
    std::vector<Particle*> particles;
private:

    //pointers to modules
    std::shared_ptr<TimeIntegration> timeIntegration;
    std::shared_ptr<DataManager> dataManager;
    std::shared_ptr<Console> console;
    std::shared_ptr<Cooling> cooling;
    std::shared_ptr<SFR> sfr;

    //calculations without the octree, just for debugging purposes
    void calculateForcesWithoutOctree(Particle* p);
};
