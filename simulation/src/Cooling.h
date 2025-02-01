#pragma once
#include <cmath>
#include <vector>
#include <memory>
#include "Particle.h"

struct CIEFractions {
        double n_H0;
        double n_Hplus;
        double n_He0;
        double n_Heplus;
        double n_Heplusplus;
        double n_e;
    };

class Cooling
{
public:
    Cooling() {}
    ~Cooling() {}

    void coolingRoutine(Particle* particle);
    
private:
    double coolingRatePrimordialCgs(
        double T,
        double n_e,    
        double n_H0,        
        double n_Hplus,    
        double n_He0,      
        double n_Heplus,   
        double n_Heplusplus,
        Particle* particle
    );

    CIEFractions getCIEFractions(
        double T,
        double n_H,
        double n_He
    );

    inline double kIonH0(double T);
    inline double kIonHe0(double T);
    inline double kIonHeplus(double T);
    inline double alphaHplus(double T);
    inline double alphaHeplus(double T);
    inline double alphaHeplusplus(double T);
    inline double alphaDielecHeplus(double T);
    inline double lambdaExcitationH0(double T);
    inline double lambdaExcitationHeplus(double T);
    inline double gauntFactorFF(double T);
};