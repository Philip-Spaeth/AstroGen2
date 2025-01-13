#include "Simulation.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include "Units.h"
#include "Log.h"

Simulation::Simulation()
{
    //construct the modules
    timeIntegration = std::make_shared<TimeIntegration>();
    dataManager = std::make_shared<DataManager>("../../output_data/");
    console = std::make_shared<Console>();
    cooling = std::make_shared<Cooling>();
    sfr = std::make_shared<SFR>();
}

Simulation::~Simulation() {}
void rotateSystemAroundX_90(std::vector<Particle*>& particles)
{
    // +90°: cos(+90°)=0, sin(+90°)=+1
    double cosA =  0.0;
    double sinA = +1.0;

    for (Particle* p : particles)
    {
        // --- Position ---
        double x = p->position.x;
        double y = p->position.y;
        double z = p->position.z;

        double yNew = y * cosA - z * sinA; // = -z
        double zNew = y * sinA + z * cosA; // = y

        p->position.x = x;
        p->position.y = yNew;
        p->position.z = zNew;

        // --- Geschwindigkeit ---
        double vx = p->velocity.x;
        double vy = p->velocity.y;
        double vz = p->velocity.z;

        double vyNew = vy * cosA - vz * sinA; // = -vz
        double vzNew = vy * sinA + vz * cosA; // = vy

        p->velocity.x = vx;
        p->velocity.y = vyNew;
        p->velocity.z = vzNew;
    }
}


bool Simulation::init()
{
    //load the config file
    if (!dataManager->loadConfig("../Config.ini", this))
    {
        std::cerr << "Error: Could not load the config file." << std::endl;
        return false;
    }
    Log::setOutputDir(dataManager->outputPath + "/logs");
    
    fixedStep = endTime / fixedTimeSteps;


//catch errors 
    //check if minTimeStep is smaller or equal to maxTimeStep
    if (minTimeStep > maxTimeStep)
    {
        std::cerr << "Error: minTimeStep is greater than maxTimeStep." << std::endl;
        return false;
    }
    //ckeck if the end time  / minTimeStep is < fixedTimeSteps
    if (endTime / minTimeStep < fixedTimeSteps)
    {
        std::cerr << "Error: endTime / minTimeStep is smaller than fixedTimeSteps." << std::endl;
        return false;
    }

    //print the computers / server computational parameters like number of threads, ram, cpu, etc.
    Console::printSystemInfo();
    
    Log::startProcess("load IC");
    dataManager->loadICs(particles, this);
    /*dataManager->inputPath = "500k_Andromeda.gal";
    std::vector<Particle*> andromedaParticles;
    dataManager->loadICs(andromedaParticles, this);
    dataManager->inputPath = "500k_Milchstraße.gal";
    std::vector<Particle*> milkyWayParticles;
    dataManager->loadICs(milkyWayParticles, this);

    double angleX = 0.0;                  // kein kippeln um X
    double angleY = 77.0 * (M_PI / 180.0);   // ~77° - Inklination
    double angleZ = 40.0 * (M_PI / 180.0);   // ~40° - Positionswinkel

    double cosX = cos(angleX), sinX = sin(angleX);
    double cosY = cos(angleY), sinY = sin(angleY);
    double cosZ = cos(angleZ), sinZ = sin(angleZ);
    for (int i = 0; i < (int)andromedaParticles.size(); i++)
    {
        vec3 pos = andromedaParticles[i]->position;

        double y1 = pos.y * cosX - pos.z * sinX;
        double z1 = pos.y * sinX + pos.z * cosX;
        double x2 = pos.x * cosY + z1 * sinY;
        double z2 = -pos.x * sinY + z1 * cosY;
        double x3 = x2 * cosZ - y1 * sinZ;
        double y3 = x2 * sinZ + y1 * cosZ;
        andromedaParticles[i]->position = vec3(x3, y3, z2);
        vec3 vel = andromedaParticles[i]->velocity;
        double vy1 = vel.y * cosX - vel.z * sinX;
        double vz1 = vel.y * sinX + vel.z * cosX;
        double vx2 = vel.x * cosY + vz1 * sinY;
        double vz2 = -vel.x * sinY + vz1 * cosY;
        double vx3 = vx2 * cosZ - vy1 * sinZ;
        double vy3 = vx2 * sinZ + vy1 * cosZ;

        andromedaParticles[i]->velocity = vec3(vx3, vy3, vz2);
        // 2.5e22m , -117000m/s, 1e22m,  time = 3091489847 years , -232643m/s
        andromedaParticles[i]->position += vec3(1e22, 0, 0);
        andromedaParticles[i]->velocity += vec3(-232643,0,0);
        particles.push_back(andromedaParticles[i]);
    }
    
    for (int i = 0; i < (int)milkyWayParticles.size(); i++)
    {
        particles.push_back(milkyWayParticles[i]);
    }

    rotateSystemAroundX_90(particles);
    dataManager->saveData(particles, 0, fixedTimeSteps, numParticlesOutput, fixedStep, endTime, 0.0);
    return false;*/
    
    //numberOfParticles = particles.size();
    if((size_t)numberOfParticles != particles.size())
    {
        std::cerr << "Error: Number of particles in the ConfigFile does not match the number of particles in the data file." << std::endl;
        std::cout << "Number of particles in the ConfigFile: " << numberOfParticles << std::endl;
        std::cout << "Number of particles in the data file: " << particles.size() << std::endl;
        return false;
    }
    
    //check if there are null pointers in the particles vector
    #pragma omp parallel for
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (!particles[i]) 
        {
            std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
        }
        if(particles[i]->type == 2)
        {
            particles[i]->T = (Constants::GAMMA - 1.0) * particles[i]->U * Constants::prtn * particles[i]->mu / (Constants::k_b);
            //std::cout << particles[i]->U << "   ,   " << particles[i]->T << std::endl;
        }
    }
    //shuffle the particles to get a random distribution
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(particles.begin(), particles.end(), g);

    //Log::saveVelocityCurve(particles, numberOfParticles);
    Log::startProcess("build tree");
    Tree* tree = new Tree(this);
    tree->buildTree();
    std::cout << "\nInitial tree size: " << std::fixed << std::scientific << std::setprecision(1) << tree->root->radius <<"m"<< std::endl;
    
    Log::startProcess("Visual Density");
    visualDensityRadius = tree->root->radius / 100000;
    //calculate the visualDensity, just for visualization
    tree->calcVisualDensity();
    //calculate the gas density for SPH
    Log::startProcess("SPH density and update");
    tree->calcGasDensity();

    // Initial force calculation
    Log::startProcess("Force Calculation");
    tree->calculateForces();
    
    //delete the tree
    Log::startProcess("delete tree");
    delete tree;
    //tree->deleteTreeParallel();

    //save the particles data
    Log::startProcess("Save data");
    dataManager->saveData(particles, 0, fixedTimeSteps, numParticlesOutput, fixedStep, endTime, 0.0);
    
    //print the memory size of the data
    double storageSize = fixedTimeSteps;
    if(dataManager->outputFormat == "ag") storageSize *= dataManager->ag_MemorySize;
    if(dataManager->outputFormat == "age") storageSize *= dataManager->age_MemorySize;
    if(dataManager->outputFormat == "agc") storageSize *= dataManager->agc_MemorySize;
    if(dataManager->outputFormat == "hdf5") storageSize *= dataManager->hdf5_MemorySize;
    if(dataManager->outputFormat == "gadget") storageSize *= dataManager->gadget_MemorySize;
    if(storageSize < 1000000000)
    {
        std::cout << std::fixed << std::setprecision(1) << "Storage size of the data: " << storageSize / 1000000 << " MB" << std::endl;
    }
    else
    {
        std::cout << std::fixed << std::setprecision(1) << "Storage size of the data: " << storageSize / 1000000000 << " GB" << std::endl;
    }

    Log::endProcess();
    return true;
}

void Simulation::run()
{
    startTimeSimulation = std::chrono::high_resolution_clock::now();


    globalTime = 0.0;
    double nextSaveTime = fixedStep;

    // Set the next integration time for each particle to 0.0 to ensure that the force is calculated in the first iteration
    #pragma omp parallel for
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (!particles[i]) 
        {
            #pragma omp critical
            std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
            continue;
        }
        particles[i]->nextIntegrationTime = 0.0;
    }

    // Initialize particles' time steps and next integration times
    #pragma omp parallel for
    for (int i = 0; i < numberOfParticles; i++)
    {
        if (!particles[i]) 
        {
            #pragma omp critical
            std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
            continue;
        }
        double accelMag = particles[i]->acc.length();
        if (accelMag > 0) {
            double timeStep = eta * std::sqrt(e0 / accelMag);
            particles[i]->timeStep = std::clamp(timeStep, minTimeStep, maxTimeStep);
            particles[i]->timeStep = std::max(std::pow(2, std::floor(std::log2(particles[i]->timeStep))), minTimeStep);
            particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
        } else {
            particles[i]->timeStep = minTimeStep;
            particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
        }
    }

    // Main simulation loop
    while (globalTime < endTime)
    {
        // Determine the next integration time for each particle
        #pragma omp parallel for
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime >= particles[i]->nextIntegrationTime)
            {
                double accelMag = particles[i]->acc.length();
                if (accelMag > 0) {
                    double timeStep = eta * std::sqrt(e0 / accelMag);
                    particles[i]->timeStep = std::clamp(timeStep, minTimeStep, maxTimeStep);
                    particles[i]->timeStep = std::max(std::pow(2, std::floor(std::log2(particles[i]->timeStep))), minTimeStep);
                    particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
                } else {
                    particles[i]->timeStep = minTimeStep;
                    particles[i]->nextIntegrationTime = globalTime + particles[i]->timeStep;
                }
            }
        }

        // Find the smallest next integration time among all particles
        double minIntegrationTime = std::numeric_limits<double>::max();
        #pragma omp parallel for reduction(min:minIntegrationTime)
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                #pragma omp critical
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (particles[i]->nextIntegrationTime < minIntegrationTime)
            {
                minIntegrationTime = particles[i]->nextIntegrationTime;
            }
        }

        // Advance global time by the smallest integration time
        globalTime = minIntegrationTime;
        
        Log::startProcess("first kick");
        // Update positions and velocities using the KDK Leapfrog scheme for particles due to be integrated
        #pragma omp parallel for
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                #pragma omp critical
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime == particles[i]->nextIntegrationTime)
            {
                timeIntegration->Kick(particles[i], particles[i]->timeStep);
                timeIntegration->Drift(particles[i], particles[i]->timeStep);
            }
        }
        
        Log::startProcess("build tree");
        Tree* tree = new Tree(this);
        tree->buildTree();

        Log::startProcess("Visual Density");
        tree->calcVisualDensity();
        Log::startProcess("SPH density and update");
        tree->calcGasDensity();
        
        Log::startProcess("Force Calculation");
        tree->calculateForces();

        // Second kick
        Log::startProcess("second kick");
        sfr->totalSFR = 0;
        #pragma omp parallel for
        for (int i = 0; i < numberOfParticles; i++)
        {
            if (!particles[i]) 
            {
                #pragma omp critical
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime == particles[i]->nextIntegrationTime)
            {
                if(particles[i]->type == 2)
                {
                    
                    //std::cout << std::fixed << std::scientific << "Particle " << i << " rho: " << particles[i]->T << std::endl;
                    //cooling and star formation
                    if(coolingEnabled)
                    {
                        cooling->coolingRoutine(particles[i]);
                    }
                    if(starFormation)
                    {
                        //calc SFR
                        sfr->sfrRoutine(particles[i]);
                    }
                    
                    if(particles[i]->type == 2)
                    {
                        // Integrate the internal energy
                        timeIntegration->Ueuler(particles[i], particles[i]->timeStep);
                    }
                }

                //aplly hubbles law
                double H0SI = (H0 * Units::KMS) / Units::MPC;
                double scale_factor = exp(H0SI * particles[i]->timeStep);
                particles[i]->position *= scale_factor;

                timeIntegration->Kick(particles[i], particles[i]->timeStep);
                // Schedule the next integration time for this particle
                particles[i]->nextIntegrationTime += particles[i]->timeStep;
            }

        }
        std::cout << "SFR: " << sfr->totalSFR << std::endl;

        Log::startProcess("delete tree");
        delete tree;
        //tree->deleteTreeParallel();

        // Save data at regular intervals defined by fixedStep
        if (globalTime >= nextSaveTime)
        {
            Log::startProcess("Save data");
            dataManager->saveData(particles, static_cast<int>(nextSaveTime / fixedStep), fixedTimeSteps, numParticlesOutput, fixedStep, endTime, globalTime);
            console->printProgress(static_cast<int>(nextSaveTime / fixedStep), fixedTimeSteps, "");
            nextSaveTime += fixedStep;
            
            if(starFormation)
            {
                double gasMass = 0;
                double totalMass = 0;
                for(int i = 0; i < numberOfParticles; i++)
                {
                    if(particles[i]->type == 2)
                    {
                        gasMass += particles[i]->mass;
                    }
                    totalMass += particles[i]->mass;
                }
                std::cout << "Gas fraction: " << gasMass / totalMass * 100 << "%" << std::endl;
            }
            if(globalTime == fixedStep * 10)
            {
                Log::avg_R_sfr(particles, numberOfParticles);
                Log::avg_R_U(particles, numberOfParticles);
            }
            Log::total_Mass(particles, globalTime);
            Log::sfr(particles, globalTime);
            Log::avg_U(particles, globalTime);
        }
    }

    std::cout << "Simulation finished." << std::endl;

    // Print the total simulation time
    auto endTimeSimulation = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTimeSimulation - startTimeSimulation);
    std::cout << "Total simulation time: " << duration.count() << " seconds" << std::endl;
}

//only for Debugging and Comparison with the octree, extremely slow
void Simulation::calculateForcesWithoutOctree(Particle* p)
{
    p->acc = vec3(0.0, 0.0, 0.0);
    p->dUdt = 0;

    #pragma omp parallel for
    for (int j = 0; j < numberOfParticles; j++)
    {
        if (p != particles[j])
        {
            vec3 d = particles[j]->position -p->position;
            double r = d.length();
            vec3 newAcceleration = d * (Constants::G * particles[j]->mass / (r * r * r));
            p->acc += newAcceleration;
        }
    }
}