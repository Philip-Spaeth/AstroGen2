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

Simulation::~Simulation(){}

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
    //dataManager->loadICs(particles, this);
    
    dataManager->inputPath = "500k_Andromeda.gal";
    dataManager->loadICs(particles, this);

    std::cout << "loaded" << std::endl;

    dataManager->saveData(particles, 0, fixedTimeSteps, numParticlesOutput, fixedStep, endTime, 0.0);

    std::cout << "saved" << std::endl;
    
    dataManager->inputPath = "../output_data/test/0.gadget_LG2";
    std::vector<Particle*> milkyWayParticles;
    dataManager->loadICs(milkyWayParticles, this);

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
    std::mt19937 g(42); 
    std::shuffle(particles.begin(), particles.end(), g);

    //Log::saveVelocityCurve(particles, numberOfParticles);
    Log::startProcess("build tree");
    Tree* tree = new Tree(this);
    tree->buildTree();
    std::cout << "\nInitial tree size: " << std::fixed << std::scientific << std::setprecision(1) << tree->root->radius <<"m"<< std::endl;
    
    visualDensityRadius = tree->root->radius / 100000;
    //calculate the visualDensity, just for visualization
    tree->calcVisualDensity();
    //calculate the gas density for SPH
    Log::startProcess("density");
    tree->calcGasDensity(N_in_h);

    // Initial force calculation
    Log::startProcess("Gravity");
    tree->calculateGravity();

    Log::startProcess("SPH Force");
    tree->calculateSPH();
    
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
    for (int i = 0; i < (int)particles.size(); i++)
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
    for (int i = 0; i < (int)particles.size(); i++)
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
        for (int i = 0; i < (int)particles.size(); i++)
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
        for (int i = 0; i < (int)particles.size(); i++)
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
        for (int i = 0; i < (int)particles.size(); i++)
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

        tree->calcVisualDensity();
        Log::startProcess("density");
        tree->calcGasDensity(N_in_h);

        Log::startProcess("SN-feedback");
        //SN pending
        if(SNFeedbackEnabled)
        {
            for (int i = 0; i < (int)particles.size(); i++)
            {
                if (!particles[i]) {continue;}
                if (particles[i]->SN_pending)
                {
                    if(particles[i]->node)
                    {
                        particles[i]->node->SNFeedback_Kawata(particles[i], 4e48, e_sn, f_v_sn);
                        particles[i]->SN_pending = false;

                        auto it = std::find(particles.begin(), particles.end(), particles[i]);
                        if (it != particles.end()) 
                        {
                            delete *it;
                            #pragma omp critical
                            {
                                particles.erase(it);
                            }
                        }
                    }
                }
            }
        }

        Log::startProcess("Gravity");
        tree->calculateGravity();

        Log::startProcess("SPH Force");
        tree->calculateSPH();

        // Second kick
        Log::startProcess("second kick");
        sfr->totalSFR = 0;
        double newStarsMass = 0;
        for (int i = 0; i < (int)particles.size(); i++)
        {
            if (!particles[i]) 
            {
                std::cerr << "Error: Particle " << i << " is not initialized." << std::endl;
                continue;
            }
            if (globalTime == particles[i]->nextIntegrationTime)
            {
                if(particles[i]->type == 2)
                {
                    //cooling and star formation
                    if(coolingEnabled)
                    {
                        cooling->coolingRoutine(particles[i]);
                    }
                    if(starFormation)
                    {
                        //calc SFR
                        sfr->sfrRoutine(particles[i], this, newStarsMass);
                        numberOfParticles = particles.size();
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

        double SFR = (newStarsMass / Units::MSUN) / (particles[23]->timeStep / Units::YR);
        sfr->totalSFR = SFR;

        Log::startProcess("delete tree");
        delete tree;

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
                for(int i = 0; i < (int)particles.size(); i++)
                {
                    if(particles[i]->type == 2)
                    {
                        gasMass += particles[i]->mass;
                    }
                    totalMass += particles[i]->mass;
                }
                std::cout << std::setprecision(6);
                std::cout << "Gas fraction: " << gasMass / totalMass * 100 << "%   Global SFR:" << sfr->totalSFR << "   N: " << particles.size() << std::endl;
                std::cout << std::fixed << std::setprecision(2);
            }
            
            if(globalTime == fixedStep * 10)
            {
                //Log::avg_R_sfr(particles, particles.size());
                //Log::avg_R_U(particles, particles.size());
            }
            Log::total_Mass(particles, globalTime);
            Log::sfr(particles, globalTime, sfr->totalSFR);
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
    for (int j = 0; j < (int)particles.size(); j++)
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