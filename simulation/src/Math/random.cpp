#include "random.h"
#include <chrono>
#include <cmath>
#include <cstdlib>


std::mt19937 random::generator(std::random_device{}());

double random::normal_dist(double mean, double stddev)
{
    std::normal_distribution<double> distribution(mean, stddev);
    return distribution(generator);

}

void random::setSeed(unsigned int seed)
{
    random::generator.seed(seed);
}

double random::uniform(double min, double max)
{
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(generator);
}