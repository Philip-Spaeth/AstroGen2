#include "kernel.h"
#include "Constants.h"

double kernel::cubicSplineKernel(double r, double h)
{
    double alpha_3d = 1.0 / (Constants::PI * h * h * h);
    double q = r / h;
    if (q < 1.0) 
    {
        return alpha_3d * (1 - 1.5 * q * q + 0.75 * q * q * q);
    }
    else if (q < 2.0) {
        return alpha_3d * 0.25 * pow(2 - q, 3);
    }
    return 0.0;
}

vec3 kernel::gradientCubicSplineKernel(const vec3& r, double h) 
{
    double r_norm = r.length();
    if (r_norm == 0.0) {
        return vec3(0.0, 0.0, 0.0); // Handle zero division gracefully
    }
    const double alpha_3d = 1.0 / (Constants::PI * h * h * h);
    double q = r_norm / h;
    double gradient_scalar;

    if (q < 1.0) {
        gradient_scalar = alpha_3d * (-3.0 * q + 2.25 * q * q);
    } else if (q < 2.0) {
        gradient_scalar = alpha_3d * (-0.75 * pow(2 - q, 2));
    } else {
        return vec3(0.0, 0.0, 0.0);
    }

    // Compute the vector gradient
    gradient_scalar /= h; // Adjust gradient by 1/h due to chain rule
    return r * (gradient_scalar / r_norm); // Normalize r and scale by gradient
}

double kernel::softeningKernel(double u)
{
    if (u >= 0 && u < 1.0/2.0)
    {
        return (16.0/3.0) *std::pow(u, 2) - (48.0/5.0) * std::pow(u, 4) + (32.0/5.0) * std::pow(u, 5) - (14.0/5.0);
    }
    else if (u >= 1.0/2.0 && u < 1.0)
    {
        return (1.0/(15.0 * u)) + (32.0/3.0) * std::pow(u, 2) - (16.0) * std::pow(u, 3) - (48.0/5.0) * std::pow(u, 4) - (32.0/15.0) * std::pow(u, 5) - (16.0/5.0);
    }
    else if (u >= 1.0)
    {
        return -1.0 / u;
    }
    return 0.0;
}