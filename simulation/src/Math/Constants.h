#pragma once

namespace Constants
{
    //Physical Constants
    constexpr double C = 299792458.0;
    constexpr double G = 6.67430e-11;

    //Mathematical Constants
    constexpr double PI = 3.14159265358979323846;
    constexpr double E = 2.71828182845904523536;
    constexpr double PHI = 1.61803398874989484820;

    //SPH related constants in SI units
    constexpr double GAMMA = 5.0/3.0; //adiabatic index for ideal gas
    constexpr double BK = 1.38064852e-23; //Boltzmann constant
    constexpr double prtn = 1.6726219e-27; //proton mass in kg
    //mean molecular weight
    constexpr double mMWH2 = 1; //neutral hydrogen
    constexpr double mMWHe = 4; //neutral helium
    constexpr double mMWHeI = 4; //ionized helium
    constexpr double mMWH2I = 0.5; //ionized hydrogen
};