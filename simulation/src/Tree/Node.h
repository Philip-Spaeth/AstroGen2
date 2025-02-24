#pragma once 

#include <iostream>
#include <memory>
#include <vector>
#include <unordered_map>
#include <future>
#include <random>
#include <cmath>
#include <queue>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdint>
#include "Particle.h"
#include "vec3.h"

class Particle;
class Node
{
private:
    //no repruduction
    bool memSafeMode = true;

public:
    Node();
    virtual ~Node();
    void deleteTreeParallel(int cores);

    void insert(const std::vector<Particle*> particles, int cores);
    int getOctant(Particle* newParticle);

    void calculateGravityForce(Particle* newparticle, double softening, double theta) const;
    void calcSPHForce(Particle* p);

    // kinematic and thermal feedback following Kawata (2001)
    void SNFeedback_Kawata(Particle* p, double snEnergy, double epsilonSN, double f_v);

    //calculate the density for all particles, without kernel approximation
    void calcVisualDensity(double radiusDensityEstimation);
    void calcDensityPart(int N, Particle* p, int type);

    //SPH only for gas particles, Stars and dark matter particles are not affected
    double gasMass = 0; //mass of gas particles
    void calcDensity(int N, Particle* p);

    int depth;
    bool isLeaf = false;
    //if leaf node, the particle is stored here
    Particle* particle;

    vec3 position;
    double radius;

    //Mass center of the node
    vec3 centerOfMass;
    double mass;
    //childparticles vector for the density calculation
    std::vector<Particle*> childParticles = std::vector<Particle*>();

    //Children of the node
    Node* children[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

    //parent of the node
    Node* parent;
};

