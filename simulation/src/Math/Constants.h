#pragma once

namespace Constants
{
    constexpr double C = 299792458.0;
    constexpr double G = 6.67430e-11;
    constexpr double PI = 3.14159265358979323846;
    constexpr double E = 2.71828182845904523536;
    constexpr double PHI = 1.61803398874989484820;

    //SPH related constants in SI units
    constexpr double GAMMA = 5.0/3.0; //adiabatic index for ideal gas
    constexpr double k_b = 1.38064852e-23; //Boltzmann constant
    constexpr double h_p = 6.62607015e-34; //Planck constant
    constexpr double prtn = 1.6726219e-27; //proton mass in kg
    const double X_H = 0.76;            // Hydrogen mass fraction
    const double X_He = 0.24;           // Helium mass fraction
};