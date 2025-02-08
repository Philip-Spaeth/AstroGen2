#pragma once 

#include <iostream>
#include "Particle.h"
#include <memory>
#include <vector>

class Particle;

class Node
{
public:
    Node();
    virtual ~Node();
    void deleteTreeParallel(int cores);

    void insert(const std::vector<Particle*> particles, int cores);
    std::vector<int> zuweiseKerne(Node* children[], size_t size, int gesamtKerne);
    //old
    void insert(Particle* newParticle);

    int getOctant(Particle* newParticle);

    void calculateGravityForce(Particle* newparticle, double softening, double theta) const;
    vec3 calcSPHForce(Particle* newparticle) const;

    void SNFeedback(Particle* p, double energy, double massInH);

    double mH = 0;
    double mRho = 0;
    double mP = 0;
    vec3 mVel = vec3(0,0,0);

    //calculate the density for all particles
    void calcVisualDensity(double radiusDensityEstimation);

    //SPH only for gas particles, Stars and dark matter particles are not affected
    double gasMass = 0; //mass of gas particles
    void calcGasDensity(double massInH);

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