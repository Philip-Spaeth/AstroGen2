#pragma once

#include "vec3.h"
#include <memory>
#include "Tree/Node.h"

class Node;

class Particle
{
public:
    Particle(vec3 position,vec3 velocity,vec3 acceleration, double mass = 1){}

    Particle(){}
    ~Particle(){}

    // Particle properties in SI units
    vec3 position = vec3(0,0,0);
    vec3 velocity = vec3(0,0,0);
    vec3 acceleration = vec3(0,0,0);
    vec3 acc = vec3(0,0,0);
    double mass = 0;

    //unique id
    unsigned int id;

    //adaptive time integration
    double timeStep = 0;
    double nextIntegrationTime = 0;

    //SPH only for gas particles, Stars and dark matter particles are not affected
    uint8_t type = 1; // 1 = star, 2 = gas, 3 = dark matter
    uint8_t galaxyPart = 1; // 1 = disk, 2 = bulge, 3 = halo

    // Fluid properties (SPH) only for gas particles
    //calculated by the Tree
    double h = 0; // smoothing length in m
    double rho = 0; //density in kg/m^3
    double P = 0; //pressure in Pascal
    double T = 0; //temperature in Kelvin
    double mu = 0.58; //mean molecular weight
    double n_H = 0; //hydrogen number density
    double n_He = 0; //helium number density
    double n_e = 0; //electron number density
    double n_H2 = 0; //molecular hydrogen number density
    double ne = 0; //free electron fraction
    double dUdt = 0;
    double U = 0; //internal energy in J/kg

    //star formation
    double sfr = 0; //star formation rate

    //calculated for all particles not just SPH gas, just for visualization
    double visualDensity = 0;

    Node* node;
};