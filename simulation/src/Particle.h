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

    double rho = 0; //total density in kg/m^3

    double rho_h = 0; //hot phase in kg/m^3
    double rho_c = 0; //cold phase in kg/m^3
    double rho_star = 0; //star phase in kg/m^3
    double rho_th = 0; // density threshold
    
    double mass_h = 0; //hot phase in kg
    double mass_c = 0; //cold phase in kg
    double mass_star = 0; //star phase in kg

    double A_evap = 0; //supernova evaporation parameter A
    double t_sfrt = 0; //star formation timescale in s

    double P = 0; //pressure in Pascal
    double T = 0; //temperature in Kelvin

    double dUdt = 0; // density threshold
    double U = 0; //internal energy in J/kg

    double u_h = 0; //internal energy hot phase in J/kg
    double du_hdt = 0; //rate of change of internal energy hot phase in J/kg/s
    double u_c = 0; //internal energy cold phase in J/kg
    double du_cdt = 0; //rate of change of internal energy cold phase in J/kg/s

    double delayedCoolingTime = 0;
    
    double mu = 0.58; //mean molecular weight
    //divergence of velocity field
    double div_v = 0;
    //sound crossing time
    double t_s = 0;

    //star formation
    double sfr = 0; //star formation rate
    bool SN_pending = false; //supernova pending

    //calculated for all particles not just SPH gas, just for visualization
    double visualDensity = 0;

    Node* node;
};