#pragma once

namespace Constants
{
    constexpr double C       = 299792458.0;             // Speed of light in vacuum [m/s]
    constexpr double G       = 6.67430e-11;              // Gravitational constant [m^3 kg^-1 s^-2]
    constexpr double PI      = 3.14159265358979323846;   // Pi
    constexpr double E       = 2.71828182845904523536;   // Euler's number
    constexpr double PHI     = 1.61803398874989484820;   // Golden ratio

    constexpr double k_b     = 1.38064852e-23;           // Boltzmann constant [J/K]
    constexpr double R       = 8.314462618;              // Ideal gas constant [J/(mol·K)]
    constexpr double GAMMA   = 5.0 / 3.0;                 // Adiabatic index (monatomic ideal gas)
    constexpr double h_p     = 6.62607015e-34;           // Planck constant [J·s]
    constexpr double hbar    = h_p / (2.0 * PI);         // Reduced Planck constant [J·s]
    constexpr double e_ch    = 1.602176634e-19;          // Elementary charge [C]
    constexpr double alpha   = 7.2973525693e-3;          // Fine-structure constant
    constexpr double prtn    = 1.6726219e-27;            // Proton mass [kg]
    constexpr double ntrn    = 1.6749275e-27;            // Neutron mass [kg]
    constexpr double elctn   = 9.10938356e-31;           // Electron mass [kg]
    constexpr double mu_0    = 1.2566370614e-6;          // Vacuum permeability [N/A²]
    constexpr double eps_0   = 8.854187817e-12;          // Vacuum permittivity [F/m]
    constexpr double NA      = 6.02214076e23;            // Avogadro constant [mol^-1]
    constexpr double sigma_T = 6.6524587321e-29;         // Thomson scattering cross-section [m²]
    constexpr double a_rad   = 7.5657e-16;               // Radiation density constant [J m^-3 K^-4]
    
    const double X_H         = 0.76;                     // Hydrogen mass fraction
    const double X_He        = 0.24;                     // Helium mass fraction

    constexpr double L_sun   = 3.828e26;                 // Solar luminosity [W]
    constexpr double R_sun   = 6.9634e8;                 // Solar radius [m]
    constexpr double T_sun   = 5772.0;                   // Solar effective temperature [K]
    constexpr double AU_m    = 1.495978707e11;           // Astronomical unit [m]
};
