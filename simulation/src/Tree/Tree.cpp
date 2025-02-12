#include "Tree.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <omp.h>

Tree::~Tree()
{
    delete root;
    root = nullptr;
}

void Tree::buildTree()
{
    //start Time
    auto start = std::chrono::high_resolution_clock::now();
    //create the root node
    root = new Node();
    //setup the root node
    root->position = vec3(0.0, 0.0, 0.0);
    root->radius = calcTreeWidth();
    root->depth = 0;
    
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);
    omp_set_nested(1);
        
    root->insert(simulation->particles, max_threads);
    
    //end Time
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    //std::cout << "Tree built in " << elapsed_seconds.count() << "s" << std::endl;
}

void Tree::calculateGravity() 
{
    const int numThreads = std::thread::hardware_concurrency();
    omp_set_num_threads(numThreads);

    int numParticles = simulation->particles.size();

   #pragma omp parallel for
    for (int i = 0; i < numParticles; ++i) 
    {
        Particle* p = simulation->particles[i];
        if(!p) continue;
        if(!root)
        {
            std::cerr << "Error: Root is not initialized." << std::endl;
            continue;
        }
        if (simulation->globalTime == p->nextIntegrationTime) 
        {
            p->acc = vec3(0.0, 0.0, 0.0);
            root->calculateGravityForce(p, simulation->e0, simulation->theta);
        }
    }
}

void Tree::calculateSPH() 
{
    const int numThreads = std::thread::hardware_concurrency();
    omp_set_num_threads(numThreads);

    int numParticles = simulation->particles.size();

    #pragma omp parallel for
    for (int i = 0; i < numParticles; ++i) 
    {
        Particle* p = simulation->particles[i];
        if(!p) continue;
        if(!root)
        {
            std::cerr << "Error: Root is not initialized." << std::endl;
            continue;
        }
        if (simulation->globalTime == p->nextIntegrationTime) 
        {
            //root->calculateSPHForce(p, simulation->e0, simulation->theta);
        }
    }

}

double Tree::calcTreeWidth()
{
    //acceptanceRatio times the standard deviation of the distances
    double acceptanceRatio = 10;

    if (simulation->numberOfParticles == 0) return 0;

    std::vector<double> distances(simulation->numberOfParticles);
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        if(!simulation->particles[i]) continue;
        distances[i] = simulation->particles[i]->position.length();
    }
    double sum = std::accumulate(distances.begin(), distances.end(), 0.0);
    double mean = sum / simulation->numberOfParticles;

    double sq_sum = std::inner_product(distances.begin(), distances.end(), distances.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / simulation->numberOfParticles - mean * mean);

    double maxAcceptableDistance = mean + acceptanceRatio * stdev;

    double max = 0;
    for (int i = 0; i < simulation->numberOfParticles; i++)
    {
        double distance = distances[i];
        if (distance <= maxAcceptableDistance && distance > max)
        {
            max = distance;
        }
    }
    return max;
}

void Tree::calcGasDensity(int N_in_h)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)simulation->particles.size(); i++)
    {
        if(simulation->particles[i] != nullptr)
        {
            if(simulation->particles[i]->type == 2)
            {
                if (simulation->particles[i]->node)
                {
                    simulation->particles[i]->node->calcDensity(N_in_h, simulation->particles[i]);
                    simulation->particles[i]->P = (Constants::GAMMA - 1.0) * simulation->particles[i]->U * simulation->particles[i]->rho;
                    simulation->particles[i]->T = (Constants::GAMMA - 1.0) * simulation->particles[i]->U * Constants::prtn * simulation->particles[i]->mu / (Constants::k_b);
                }
            }
        }
    }
}

void Tree::calcVisualDensity()
{
    const int numParticles = simulation->numberOfParticles;
    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++)
    { 
        if(!simulation->particles[i]) continue;
        simulation->particles[i]->visualDensity = 0;
    }

    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++)
    { 
        if(!simulation->particles[i]) continue;
        if(simulation->particles[i]->visualDensity != 0)
        {
            continue;
        }
        if (simulation->particles[i]->node)
        {
            simulation->particles[i]->node->calcVisualDensity(simulation->visualDensityRadius);
        }
    }
}