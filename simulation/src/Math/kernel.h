#pragma once
#include "vec3.h"

class kernel
{
public:
	//SPH Kernels
	static double cubicSplineKernel(double r, double h);
	static vec3 gradientCubicSplineKernel(const vec3& r, double h);

	//softening kernel described by Springel, Yoshida & White (2001) eq. 71
	static double softeningKernel(double u);
};