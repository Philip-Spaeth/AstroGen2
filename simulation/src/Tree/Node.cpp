#include "Node.h"
#include "Constants.h"
#include "kernel.h"
#include <algorithm>
#include "omp.h"
#include <thread>
#include <iostream>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <future>
#include <random>
#include <cmath>
#include <queue>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdint>

Node::Node()
{
    mass = 0.0;
    centerOfMass = vec3(0.0, 0.0, 0.0);
    isLeaf = true;
    particle = nullptr;
    for (int i = 0; i < 8; i++)
    {
        if (children[i] != nullptr) {
            delete children[i];
            children[i] = nullptr;
        }
    }
    parent = nullptr;

    radius = 0;
}

Node::~Node()
{
    //#pragma omp parallel for
    for (int i = 0; i < 8; i++)
    {

        if (children[i])
        {
            if(memSafeMode)
            {
                if(reinterpret_cast<std::uintptr_t>(children[i]) < 0x1000) {
                    std::cerr << "Error (Dekonstruktor): children[i] pointer is invalid (address: " << children[i] << ")" << std::endl;
                    children[i] = nullptr;
                    continue;
                }
            }


            delete children[i];
            children[i] = nullptr;
        }
    }
}

void Node::deleteTreeParallel(int cores)
{
    if(cores > 1)
    {
        #pragma omp parallel for num_threads(cores)
        for (int i = 0; i < 8; i++)
        {
            if (children[i])
            {
                if(memSafeMode)
                {
                    if(reinterpret_cast<std::uintptr_t>(children[i]) < 0x1000) {
                        std::cerr << "Error (Dekonstruktor): children[i] pointer is invalid (address: " << children[i] << ")" << std::endl;
                        children[i] = nullptr;
                        continue;
                    }
                }

                children[i]->deleteTreeParallel((cores / 8) - 1);
                children[i] = nullptr;
            }
        }
    }
    else
    {
        for (int i = 0; i < 8; i++)
        {
            if (children[i])
            {
                delete children[i];
                children[i] = nullptr;
            }
        }
    }
}

// kinematic and thermal feedback following Kawata (2001)
void Node::SNFeedback_Kawata(Particle* p, double snEnergy, double epsilonSN, double f_v)
{
    if(p->h == 0) return;
    if(parent != nullptr)
    {
        if(memSafeMode)
        {
            if(reinterpret_cast<std::uintptr_t>(parent) < 0x100000) 
            {
                std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
                return;
            }
        }
        if(radius < p->h * 5)
        {
            parent->SNFeedback_Kawata(p, snEnergy, epsilonSN, f_v);
        }
    }

    if(radius >= p->h * 5)
    {
        double E_SN_i = snEnergy * epsilonSN;

        for (int i = 0; i < (int)childParticles.size(); i++)
        {
            if(childParticles[i]->type == 2)
            {
                vec3 d = childParticles[i]->position - p->position;
                double r = d.length();
                if(r < p->h)
                {
                    double w = kernel::cubicSplineKernel(r, p->h);
                    // \Delta E_{SN,j} = E_{SN,i} * (m_j / rho_g,i) * W(r_{ij}, h_{ij})
                    //    star.gasDensity ~ rho_{g,i}
                    double deltaESN_j = E_SN_i * (childParticles[i]->mass / childParticles[i]->rho) * w;

                    //kineitc fraction
                    double deltaESN_kin = f_v * deltaESN_j;
                    double deltaV_mag = std::sqrt( 2.0 * deltaESN_kin / childParticles[i]->mass);
                    //std::cout << "deltaV_mag: " << deltaV_mag << "    j: " << deltaV_mag / childParticles[i]->velocity.length() << std::endl;
                    vec3 dir = (r > 0.0) ? (d / r) : vec3{0.0, 0.0, 0.0};
                    vec3 deltaV = deltaV_mag * dir;
                    childParticles[i]->velocity += deltaV;

                    //thermal fraction
                    double deltaESN_therm = (1.0 - f_v) * deltaESN_j;

                    // \Delta U_{j} = \Delta E_{SN,j} / m_j
                    double deltaU = deltaESN_therm / childParticles[i]->mass;
                    //std::cout << "deltaU: " << deltaU << " & :  " << deltaU / childParticles[i]->U << std::endl;
                    childParticles[i]->U += deltaU;
                }
            }
        }
    }
}

void Node::calcSPHForce(Particle* p)
{
    if(!p) return;
    if(p->h == 0) return;
    if(radius == 0) return;
    if(p->type != 2) return;
    if(radius >  2* p->h)
    {
        for(int i = 0; i < (int)childParticles.size(); i++)
        {
            if(childParticles[i]->type == 2)
            {
                vec3 d = childParticles[i]->position - p->position;
                double r = d.length();
                if(r < p->h)
                {
                    vec3 acc = vec3(0,0,0);

                    vec3 v_i = p->velocity;
                    vec3 v_j = childParticles[i]->velocity;
                    vec3 v_ij = v_i - v_j;

                    double h_i = p->h;
                    double h_j = childParticles[i]->h;
                    double h_ij = (h_i + h_j) / 2.0;

                    double rho_i = p->rho;
                    double rho_j = childParticles[i]->rho;

                    double P_i = p->P;
                    double P_j = childParticles[i]->P;

                    double c_i = sqrt(Constants::GAMMA * P_i / rho_i);
                    double c_j = sqrt(Constants::GAMMA * P_i / rho_i);
                    double c_ij = (c_i + c_j) / 2.0;

                    //Pressure force
                    //Monaghan (1992)
                    acc += - childParticles[i]->mass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j)) * kernel::gradientCubicSplineKernel(d, h_i);

                    //Artificial viscosity
                    //Monaghan & Gingold (1983)
                    double MU_ij = 0.0;
                    double alpha = 0.5;
                    double beta = 1;
                    double eta = 0.01;
                    double mu_ij = h_ij * v_ij.dot(d) / (r * r + eta * (h_ij * h_ij));
                    if(v_ij.dot(d) < 0)
                    {
                        MU_ij = -alpha * c_ij * mu_ij + beta * (mu_ij * mu_ij);
                    }
                    acc += -childParticles[i]->mass * MU_ij * kernel::gradientCubicSplineKernel(d, h_ij);

                    //Internal energy
                    p->dUdt += 1.0 / 2.0 * childParticles[i]->mass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j) + MU_ij) * v_ij.dot(kernel::gradientCubicSplineKernel(d, h_i));
                    //std::cout << radius << "  ,  " << p->dUdt << " , " << childParticles.size() << std::endl;

                    if(std::isnan(acc.x) || std::isnan(acc.y) || std::isnan(acc.z)) acc = vec3(0,0,0);

                    p->acc += acc;
                }
            }
        }
    }
    else
    {
        if(parent == nullptr) return;
        if(memSafeMode)
        {
            if(reinterpret_cast<std::uintptr_t>(parent) < 0x100000) 
            {
                std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
                return;
            }
        }
        parent->calcSPHForce(p);
    }
}

void Node::calculateGravityForce(Particle* newparticle, double softening, double theta) const
{
    if (mass == 0) {
        return;
    }
    if (newparticle == nullptr) {
        return;
    }
    if (newparticle == this->particle) {
        return;
    }
    if (newparticle->mass == 0) {
        return;
    }
    vec3 d = centerOfMass - newparticle->position;
    double r = d.length();

    if (r == 0) {
        return;
    }
    
    if (isLeaf)
    {
        if (this->particle != nullptr && newparticle != this->particle)
        {
            double e0 = softening;
            
            double softeningFactorInput = r / (2.8 * e0);

            double kernelValue = kernel::softeningKernel(softeningFactorInput);
            if ((kernelValue - r) == 0) {
                return;
            }

            double e = -(2.8 * e0) / (kernelValue - r);
            if(std::isnan(e))
            {
                std::cout << "NAN e" << std::endl;
                return;
            }
            vec3 gravityAcceleration = (Constants::G * mass / (r * r + e0 * e0)) * d.normalize();
            if(abs(e) > (e0 * 0.0001) && abs(e) < 1e30)
            {
                gravityAcceleration = (Constants::G * mass / (r * r + e * e)) * d.normalize();
            }
            
            if(std::isnan(gravityAcceleration.x) || std::isnan(gravityAcceleration.y) || std::isnan(gravityAcceleration.z))
            {
                std::cout << "NAN gravacc" << std::endl;
                return;
            }
            
            newparticle->acc += gravityAcceleration;
        }
    }
    else
    {
        double s = radius / r;

        if (s < theta)
        {
            double e0 = softening;
            
            double softeningFactorInput = r / (2.8 * e0);

            double kernelValue = kernel::softeningKernel(softeningFactorInput);
            if ((kernelValue - r) == 0) {
                return;
            }

            double e = -(2.8 * e0) / (kernelValue - r);
            if(std::isnan(e))
            {
                std::cout << "NAN e" << std::endl;
                return;
            }
            vec3 gravityAcceleration = (Constants::G * mass / (r * r + e0 * e0)) * d.normalize();
            if(abs(e) > (e0 * 0.0001) && abs(e) < 1e30)
            {
                gravityAcceleration = (Constants::G * mass / (r * r + e * e)) * d.normalize();
            }
            
            if(std::isnan(gravityAcceleration.x) || std::isnan(gravityAcceleration.y) || std::isnan(gravityAcceleration.z))
            {
                std::cout << "NAN gravacc" << std::endl;
                return;
            }
            
            newparticle->acc += gravityAcceleration;
        }
        else
        { 
            for (int i = 0; i < 8; i++)
            {
                if (children[i] == nullptr) {
                    continue;
                }

                if (children[i]->mass == 0) {
                    continue;
                }

                children[i]->calculateGravityForce(newparticle, softening, theta);
            }
        }
    }
}

void Node::insert(const std::vector<Particle*> particles, int cores)
{
    if(particles.size() == 0) return;

    if(particles.size() == 1)
    {
        isLeaf = true;
        this->particle = particles[0];
        this->particle->node = this;
        centerOfMass = particle->position;
        mass = particle->mass;
        gasMass = (particle->type == 2) ? particle->mass : 0.0;
        return;
    }

    isLeaf = false;

    for (int i = 0; i < 8; i++)
    {
        children[i] = new Node();
        children[i]->position = position + vec3(
            radius * (i & 1 ? 0.5 : -0.5),
            radius * (i & 2 ? 0.5 : -0.5),
            radius * (i & 4 ? 0.5 : -0.5));
        children[i]->radius = radius / 2;
        children[i]->depth = depth + 1;
        children[i]->parent = this;
    }

    double total_mass = 0.0;
    double total_gasMass = 0.0;
    double total_position_mass_x = 0.0;
    double total_position_mass_y = 0.0;
    double total_position_mass_z = 0.0;

    int max_threads = cores;
    std::vector<std::vector<std::vector<Particle*>>> thread_octants(max_threads, std::vector<std::vector<Particle*>>(8));

    int chunk_size = static_cast<int>(particles.size() / (max_threads * (depth + 1)));
    chunk_size = std::max(chunk_size, 1);

    #pragma omp parallel num_threads(max_threads) default(none) shared(particles, thread_octants, chunk_size, max_threads) \
        reduction(+:total_mass, total_gasMass, total_position_mass_x, total_position_mass_y, total_position_mass_z)
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for schedule(dynamic, chunk_size)
        for (size_t i = 0; i < particles.size(); ++i) {

            if(depth == 0) particles[i]->node = nullptr;

            const auto& p = particles[i];
            if (!p) continue;

            total_mass += p->mass;
            if (p->type == 2) {
                total_gasMass += p->mass;
                #pragma omp atomic
                gasMass += p->mass;
            }


            total_position_mass_x += p->position.x * p->mass;
            total_position_mass_y += p->position.y * p->mass;
            total_position_mass_z += p->position.z * p->mass;

            int octant = getOctant(p);
            if (octant != -1 && thread_id < max_threads) {
                thread_octants[thread_id][octant].push_back(p);
            }
        }
    }

    mass = total_mass;
    gasMass = total_gasMass;
    if (mass > 0.0) {
        centerOfMass = vec3{ total_position_mass_x / mass, total_position_mass_y / mass, total_position_mass_z / mass };
    }

    #pragma omp parallel for schedule(static) default(none) shared(thread_octants, children, max_threads)
    for (int o = 0; o < 8; ++o) {
        for (int t = 0; t < max_threads; ++t) {
            children[o]->childParticles.insert(children[o]->childParticles.end(),
                thread_octants[t][o].begin(),
                thread_octants[t][o].end());
        }
    }


    for (int i = 0; i < 8; i++)
    {
        if (children[i]->childParticles.size() > 0)
        {
            children[i]->insert(children[i]->childParticles, cores);
        }
    }

}

int Node::getOctant(Particle* newParticle) 
{
    if(!newParticle) return -1;

    if (newParticle->position.x < position.x - radius || newParticle->position.x > position.x + radius ||
        newParticle->position.y < position.y - radius || newParticle->position.y > position.y + radius ||
        newParticle->position.z < position.z - radius || newParticle->position.z > position.z + radius) {
        return -1;
    }

    int octant = 0;
    if (newParticle->position.x > position.x) octant |= 1;
    if (newParticle->position.y > position.y) octant |= 2;
    if (newParticle->position.z > position.z) octant |= 4;
    
    return octant;
}

void Node::calcDensity(int N, Particle* p)
{
    if (childParticles.size() <= (size_t)N)
    {
        if (parent != nullptr)
        {
            if (memSafeMode && reinterpret_cast<std::uintptr_t>(parent) < 0x100000) {
                std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
                return;
            }
            parent->calcDensity(N, p);
        }
        return;
    }

    int gasN = 0;
    for (auto* particle : childParticles)
    {
        gasN += (particle->type == 2);
    }

    if (gasN < N)
    {
        if (parent != nullptr)
        {
            if (memSafeMode && reinterpret_cast<std::uintptr_t>(parent) < 0x100000) {
                std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
                return;
            }
            parent->calcDensity(N, p);
        }
        return;
    }

    using ParticleDist = std::pair<double, Particle*>;
    std::vector<ParticleDist> distances;

    for (auto* particle : childParticles)
    {
        if (particle == p || particle->type != 2) continue;
        distances.emplace_back((particle->position - p->position).length(), particle);
    }

    if ((int)distances.size() > N)
    {
        std::nth_element(distances.begin(), distances.begin() + N, distances.end());
        distances.resize(N);
    }

    double maxDistance = 0;
    for (const auto& [dist, particle] : distances)
    {
        maxDistance = std::max(maxDistance, dist);
    }
    p->h = maxDistance;

    p->rho = 0;
    for (const auto& [dist, particle] : distances)
    {
        p->rho += particle->mass * kernel::cubicSplineKernel(dist, p->h);
    }
}

void Node::calcVisualDensity(double radiusDensityEstimation) 
{
    if(parent != nullptr) {
        if(memSafeMode)
        {
        if(reinterpret_cast<std::uintptr_t>(parent) < 0x100000) {
            std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
            return;
        }
        if(reinterpret_cast<std::uintptr_t>(parent) == 0x4433c6b8197e3a36 || reinterpret_cast<std::uintptr_t>(parent) == 0x445a2118c02e3386) {
            std::cerr << "Error: specific adress (address: " << parent << ")" << std::endl;
            return;
        }
        }   

        double radiusDifference = radiusDensityEstimation - radius;
        double parentRadius = parent->radius;
        double parentRadiusDifference = radiusDensityEstimation - parentRadius;

        if(std::abs(radiusDifference) > std::abs(parentRadiusDifference) && parent) {
            parent->calcVisualDensity(radiusDensityEstimation);
        }
        else {
            double volume = radius * radius * radius;
            if(volume == 0 || mass == 0) return;
            double density = mass / volume;

            if(density == 0 || density == INFINITY) return;

            size_t numChildren = childParticles.size();
            for (size_t i = 0; i < numChildren; ++i) {
                Particle* currentParticle = childParticles[i];
                if(memSafeMode)
                { 
                if (reinterpret_cast<std::uintptr_t>(currentParticle) < 0x100000) {
                    // Handle invalid particle address
                    continue;
                }
                }

                childParticles[i]->visualDensity = density;
            }
        }
    }
}
