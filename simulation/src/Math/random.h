#pragma once

#include "vec3.h"
#include <random>


class random
{
public:
	//random Functions
	static void setSeed(unsigned int seed);
	static double uniform(double min, double max);
	static double normal_dist(double mean = 0.0, double stddev = 1.0);
private:
    static std::mt19937 generator;
};