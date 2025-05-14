#pragma once

#pragma once

namespace Units
{
    // --- Length Units (in meters) ---
    constexpr double M     = 1.0;                      // meter
    constexpr double KM    = 1.0e3;                    // kilometer
    constexpr double AU    = 1.495978707e11;           // astronomical unit
    constexpr double LY    = 9.4607e15;                // light-year
    constexpr double PC    = 3.08567758149137e16;      // parsec
    constexpr double KPC   = 1.0e3 * PC;               // kiloparsec
    constexpr double MPC   = 1.0e6 * PC;               // megaparsec
    constexpr double GPC   = 1.0e9 * PC;               // gigaparsec

    // --- Time Units (in seconds) ---
    constexpr double S     = 1.0;                      // second
    constexpr double MIN   = 60.0;                     // minute
    constexpr double H     = 3600.0;                   // hour
    constexpr double DAY   = 86400.0;                  // day
    constexpr double YR    = 3.15576e7;                // Julian year
    constexpr double KYR   = 1.0e3 * YR;               // kiloyear
    constexpr double MYR   = 1.0e6 * YR;               // megayear
    constexpr double GYR   = 1.0e9 * YR;               // gigayear

    // --- Mass Units (in kilograms) ---
    constexpr double KG     = 1.0;                     // kilogram
    constexpr double G      = 1.0e-3;                  // gram
    constexpr double MG     = 1.0e-6;                  // milligram
    constexpr double MSUN   = 1.98847e30;              // solar mass
    constexpr double MEARTH = 5.9722e24;               // Earth mass
    constexpr double MJUP   = 1.89813e27;              // Jupiter mass

    // --- Velocity Units (in meters per second) ---
    constexpr double MPS    = 1.0;                     // meter per second
    constexpr double KMS    = 1.0e3;                   // kilometer per second

    // --- Energy Units (in joules) ---
    constexpr double J      = 1.0;                     // joule
    constexpr double ERG    = 1.0e-7;                  // erg
    constexpr double EV     = 1.60218e-19;             // electronvolt
    constexpr double KEV    = 1.0e3 * EV;              // kiloelectronvolt
    constexpr double MEV    = 1.0e6 * EV;              // megaelectronvolt
    
    constexpr double N      = 1.0;                     // newton
    constexpr double PA     = 1.0;                     // pascal
    constexpr double BAR    = 1.0e5;                   // bar
    constexpr double CM3    = 1.0e-6;                  // cubic centimeter in m^3
    constexpr double M3     = 1.0;                     // cubic meter
    constexpr double GCM3   = 1.0e3;                   // g/cm³ in kg/m³
    constexpr double CM2    = 1.0e-4;                  // square centimeter
    constexpr double M2     = 1.0;                     // square meter
    constexpr double K      = 1.0;                     // kelvin
    constexpr double L      = 1.0e-3;                  // liter in m³
    constexpr double c_conv = 7.297738399394e-52;
};