#include "SFR.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "Constants.h"
#include <random>
#include "Simulation.h"
#include <algorithm> 

SFR::SFR(){rng.seed(42);}

void SFR::sfrRoutine(Particle* particle, Simulation* sim, double& newStarMass)
{
    if(!particle) return;
    if(particle->rho <= 0) return;
    if(particle->U <= 0) return;
    if(particle->mass <= 0) return;

    return;

    double dt = particle->timeStep;
    
    double t_sfr0 = 8e16;

    double beta = 0.1;
    double A0 = 1000;

    double u_c = particle->u_c;
    double u_sn = 1e51;
    double u_4 = 4e48;

    //lamba of (rho, u_sn/A0)
    double lamba_net = 1e-25;
    double lambda = lamba_net / pow(particle->rho, 2);

    double x_th = 1 - (A0 * (u_4 / u_sn));
    double rho_th = (x_th / pow((1-x_th), 2)) * (beta * u_sn - ((1-beta) * u_c)) / (t_sfr0 * lambda);

    //calc t_sfrt
    double t_sfr = t_sfr0 * pow(particle->rho / rho_th, -0.5);

    double m_sf = particle->mass_c * (dt / particle->t_sfrt);
}