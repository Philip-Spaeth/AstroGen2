#include "Node.h"
#include "Constants.h"
#include "kernel.h"
#include <algorithm>
#include "omp.h"
#include <thread>
#include <iostream>
#include <numeric>
#include <vector>
#include <future>

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
            if(reinterpret_cast<std::uintptr_t>(children[i]) < 0x1000) {
                std::cerr << "Error (Dekonstruktor): children[i] pointer is invalid (address: " << children[i] << ")" << std::endl;
                children[i] = nullptr;
                continue;
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
                if(reinterpret_cast<std::uintptr_t>(children[i]) < 0x1000) {
                    std::cerr << "Error (Dekonstruktor): children[i] pointer is invalid (address: " << children[i] << ")" << std::endl;
                    children[i] = nullptr;
                    continue;
                }

                children[i]->deleteTreeParallel((cores / 8) - 1);
                //delete children[i];
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



vec3 Node::calcSPHForce(Particle* newparticle) const
{
    vec3 acc = vec3(0,0,0);
    double h_i = newparticle->h;
    double h_j = mH;
    if(isLeaf) h_j = this->particle->h;
    h_j = h_i;
    double h_ij = (h_i + h_j) / 2.0;
    //if(h_i == 0 || h_j == 0) return vec3(0,0,0);

    double rho_i = newparticle->rho;
    double rho_j = mRho;
    if(isLeaf) rho_j = this->particle->rho;
    rho_j = rho_i;
    double rho_ij = (rho_i + rho_j) / 2.0;
    //if(rho_i == 0 || rho_j == 0) return vec3(0,0,0);

    double P_i = newparticle->P;
    double P_j = mP;
    if(isLeaf) P_j = this->particle->P;
    P_j = P_i;
    //if(P_i == 0 || P_j == 0) return vec3(0,0,0);

    vec3 v_i = newparticle->velocity;
    vec3 v_j = mVel;
    if(isLeaf) v_j = this->particle->velocity;
    vec3 v_ij = v_i - v_j;

    vec3 d = newparticle->position - centerOfMass;
    double r = d.length();
    double c_i = sqrt(Constants::GAMMA * P_i / rho_i);
    double c_j = sqrt(Constants::GAMMA * P_i / rho_i);
    double c_ij = (c_i + c_j) / 2.0;

    //Pressure force
    //Monaghan (1992)
    //std::cout << d << std::endl;
    acc += -gasMass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j)) * kernel::gradientCubicSplineKernel(d, h_i);
    //Artificial viscosity
    //Monaghan & Gingold (1983)
    double MU_ij = 0.0;
    double alpha = 0.2;
    double beta = 0.4;
    double eta = 0.01;
    double mu_ij = h_ij * v_ij.dot(d) / (r * r + eta * (h_ij * h_ij));
    if(v_ij.dot(d) < 0)
    {
        MU_ij = -alpha * c_ij * mu_ij + beta * (mu_ij * mu_ij);
    }
    //std::cout << -gasMass * MU_ij * kernel::gradientCubicSplineKernel(d, h_ij) << std::endl;
    acc += -gasMass * MU_ij * kernel::gradientCubicSplineKernel(d, h_ij);

    //Internal energy
    newparticle->dUdt += 1.0 / 2.0 * gasMass * (P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j) + MU_ij) * v_ij.dot(kernel::gradientCubicSplineKernel(d, h_i));

    if(std::isnan(acc.x) || std::isnan(acc.y) || std::isnan(acc.z)) return vec3(0,0,0);

    return acc;
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

            if (r < newparticle->h * 2)
            {
                if (this->particle->type == 2 && newparticle->type == 2)
                {
                    vec3 sphForce = calcSPHForce(newparticle);
                    newparticle->acc += sphForce;
                }
            }
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

            if (r < newparticle->h * 2)
            {
                if (newparticle->type == 2 && gasMass > 0)
                {
                    vec3 sphForce = calcSPHForce(newparticle);
                    newparticle->acc += sphForce;
                }
            }
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

    if((int)particles.size() < cores*100)
    {
        for (size_t i = 0; i < particles.size(); i++)
        {
            if (!particles[i]) continue;
            insert(particles[i]);
        }
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

    int max_threads = cores > 0 ? cores : omp_get_max_threads();
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
                gasMass += p->mass;
                if (gasMass > 0) 
                {
                    mVel = (mVel * (gasMass - p->mass) + p->velocity * p->mass) / gasMass;
                    /*if (newParticle->rho > 0) {
                        mRho = (mRho * (gasMass - newParticle->mass) + newParticle->rho * newParticle->mass) / gasMass;
                    }
                    if (newParticle->P > 0) {
                        mP = (mP * (gasMass - newParticle->mass) + newParticle->P * newParticle->mass) / gasMass;
                    }
                    if (newParticle->h > 0) {
                        mH = (mH * (gasMass - newParticle->mass) + newParticle->h * newParticle->mass) / gasMass;
                    }*/
                }
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



std::vector<int> Node::zuweiseKerne(Node* children[], size_t size, int gesamtKerne) {
    long long gesamteArbeitslast = 0;

    #pragma omp parallel for reduction(+:gesamteArbeitslast)
    for (size_t i = 0; i < size; ++i) {
        if (children[i]) {
            gesamteArbeitslast += children[i]->childParticles.size();
        }
    }

    if (gesamteArbeitslast == 0) {
        return std::vector<int>(size, 0);
    }

    std::vector<int> zuweisungen(size, 0);
    std::vector<std::pair<double, int>> gebrocheneTeile(size);

    #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        if (children[i]) {
            double exakt = (static_cast<double>(children[i]->childParticles.size()) / gesamteArbeitslast) * gesamtKerne;
            zuweisungen[i] = static_cast<int>(std::floor(exakt));
            double gebrochen = exakt - std::floor(exakt);
            gebrocheneTeile[i] = {gebrochen, static_cast<int>(i)};
        } else {
            zuweisungen[i] = 0;
            gebrocheneTeile[i] = {0.0, static_cast<int>(i)};
        }
    }

    int sumZuweisungen = std::accumulate(zuweisungen.begin(), zuweisungen.end(), 0);
    int verbleibend = gesamtKerne - sumZuweisungen;

    if (verbleibend > 0) {
        std::sort(gebrocheneTeile.begin(), gebrocheneTeile.end(),
                  [](const std::pair<double, int>& a, const std::pair<double, int>& b) -> bool {
                      return a.first > b.first;
                  });

        for (const auto& teil : gebrocheneTeile) {
            if (verbleibend <= 0)
                break;
            zuweisungen[teil.second] += 1;
            verbleibend -= 1;
        }
    }

    return zuweisungen;
}


void Node::insert(Particle* newParticle) 
{
    if (!newParticle) 
    {
        std::cerr << "Error: Particle to insert is null" << std::endl;
        return;
    }

    // Check if the particle is within the bounds of the current node
    if (newParticle->position.x < position.x - radius || newParticle->position.x > position.x + radius ||
        newParticle->position.y < position.y - radius || newParticle->position.y > position.y + radius ||
        newParticle->position.z < position.z - radius || newParticle->position.z > position.z + radius) 
    {
        //std::cerr << "Warning: Particles are outside the bounds and will not be properly considered" << std::endl;
        return;
    }
    
    //add the particle to the childParticles vector
    childParticles.push_back(newParticle);

    // if the node is a leaf node
    if (isLeaf) 
    {
        // if the node has no particle
        if (!particle) {
            // set the particle to the node
            particle = newParticle;
            particle->node = this;
        }
        // if the node has a particle
        else 
        {
            // create the children of the node
            for (int i = 0; i < 8; i++) 
            {
                children[i] = new Node();
                // setup the node properties
                children[i]->position = position + vec3(
                    radius * (i & 1 ? 0.5 : -0.5),
                    radius * (i & 2 ? 0.5 : -0.5),
                    radius * (i & 4 ? 0.5 : -0.5));
                // std::cout << children[i]->position << std::endl;
                children[i]->radius = radius / 2;
                children[i]->depth = depth + 1;
                children[i]->parent = this;
            }

            // get the octant of the particle
            int octant = getOctant(particle);
            
            // insert the particle in the corresponding octant
            if (octant != -1) 
            {
                children[octant]->insert(particle);
            }
            // get the octant of the new particle
            // insert the particle in the corresponding octant
            // get the octant of the particle
            octant = getOctant(newParticle);
            if (octant != -1) 
            {
                children[octant]->insert(newParticle);
            }

            // set the node as not a leaf node
            isLeaf = false;
            // set the particle to null
            particle = nullptr;
        }
    }
    // if the node is not a leaf node
    else 
    {
        // get the octant of the particle
        int octant = getOctant(newParticle);
        if (octant != -1) 
        {
            children[octant]->insert(newParticle);
        }
    }
    // Update the center of mass and mass
    mass += newParticle->mass;
    if(newParticle->type == 2)
    {
        gasMass += newParticle->mass;
        //calc m Node properties
        if (gasMass > 0) 
        {
            mVel = (mVel * (gasMass - newParticle->mass) + newParticle->velocity * newParticle->mass) / gasMass;
            /*if (newParticle->rho > 0) {
                mRho = (mRho * (gasMass - newParticle->mass) + newParticle->rho * newParticle->mass) / gasMass;
            }
            if (newParticle->P > 0) {
                mP = (mP * (gasMass - newParticle->mass) + newParticle->P * newParticle->mass) / gasMass;
            }
            if (newParticle->h > 0) {
                mH = (mH * (gasMass - newParticle->mass) + newParticle->h * newParticle->mass) / gasMass;
            }*/
        }
    }

    centerOfMass = (centerOfMass * (mass - newParticle->mass) + newParticle->position * newParticle->mass) / mass;
}


int Node::getOctant(Particle* newParticle) 
{
    if(!newParticle) return -1;

    if (newParticle->position.x < position.x - radius || newParticle->position.x > position.x + radius ||
        newParticle->position.y < position.y - radius || newParticle->position.y > position.y + radius ||
        newParticle->position.z < position.z - radius || newParticle->position.z > position.z + radius) {
        //std::cout << "Particle is outside the bounds" << std::endl;
        return -1; // Particle is outside the bounds
    }

    int octant = 0;
    if (newParticle->position.x > position.x) octant |= 1;
    if (newParticle->position.y > position.y) octant |= 2;
    if (newParticle->position.z > position.z) octant |= 4;
    
    return octant;
}

//consider only the gas particles and th gasMass
void Node::calcGasDensity(double massInH)
{
    if(gasMass == 0) return;

    //check if the parent ist not expired
    if(parent != nullptr)
    {

        if(reinterpret_cast<std::uintptr_t>(parent) < 0x100000) { // Beispiel für eine ungültige Adresse
            std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
            // Optional: Fehlerbehandlung, z.B. Rückkehr oder Ausnahme werfen
            return;
        }

        //go upwards the tree until the mass of the node is closest to  massInH
        if(gasMass <  massInH && parent)	//if the mass of the node is smaller than massInH and the node has a parent
        {
            //stop if the current gasMass is closer to massInH than the parent gasMass
            double massDifference =  massInH - gasMass;
            double gasMass = parent->gasMass;
            double parentMassDifference =  massInH - gasMass;
            if(std::abs(massDifference) > std::abs(parentMassDifference))
            {
                parent->calcGasDensity( massInH);
            }
        }
        //stop if the current gasMass is closer to  massInH than the parent gasMass
        double massDifference =  massInH - gasMass;
        double parentMassDifference =  massInH - parent->gasMass;
        if(std::abs(massDifference) < std::abs(parentMassDifference) && gasMass != 0)
        {
            //calculate h for all the child particles in the node
            #pragma omp parallel for 
            for (size_t i = 0; i < childParticles.size(); i++)
            {
                if(reinterpret_cast<std::uintptr_t>(childParticles[i]) < 0x1000) { // Beispiel für eine ungültige Adresse
                    std::cerr << "Error (calcGasDensity): childParticle pointer is invalid (address: " << childParticles[i] << ")" << std::endl;
                    // Optional: Fehlerbehandlung, z.B. Rückkehr oder Ausnahme werfen
                    continue;
                }

                if(childParticles[i]->type == 2)
                {
                    childParticles[i]->h = radius * 2;
                }
            }

            //calculate the rho for all the particles in the node
            double rho = 0;
            #pragma omp parallel for
            for(size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    double drho = childParticles[i]->mass * kernel::cubicSplineKernel(vec3(childParticles[i]->position - centerOfMass).length(), childParticles[i]->h);
                    rho += drho;
                }
            }

            //add the new rho to all the particles in the node
            #pragma omp parallel for
            for(size_t i = 0; i < childParticles.size(); i++)
            {
                if(childParticles[i]->type == 2)
                {
                    childParticles[i]->rho = rho;
                    //calc P, P = (gamma-1)*u*rho
                    childParticles[i]->P = (Constants::GAMMA - 1.0) * childParticles[i]->U * childParticles[i]->rho;
                    //calc T, T = (gamma-1)*u*prtn / (bk)
                    childParticles[i]->T = (Constants::GAMMA - 1.0) * childParticles[i]->U * Constants::prtn * childParticles[i]->mu / (Constants::k_b);
                }
            }
        }
    }
}

void Node::calcVisualDensity(double radiusDensityEstimation) 
{
    if(parent != nullptr) {

        // Additional check for the validity of the parent pointer
        if(reinterpret_cast<std::uintptr_t>(parent) < 0x100000) { // Example invalid address
            std::cerr << "Error: parent pointer is invalid (address: " << parent << ")" << std::endl;
            return;
        }
        if(reinterpret_cast<std::uintptr_t>(parent) == 0x4433c6b8197e3a36 || reinterpret_cast<std::uintptr_t>(parent) == 0x445a2118c02e3386) { // Example invalid address
            std::cerr << "Error: specific adress (address: " << parent << ")" << std::endl;
            return;
        }


        // Potential Issue: 'radius' is used before it's defined
        double radiusDifference = radiusDensityEstimation - radius;
        double parentRadius = parent->radius; // 'radius' is defined here
        double parentRadiusDifference = radiusDensityEstimation - parentRadius;

        if(std::abs(radiusDifference) > std::abs(parentRadiusDifference) && parent) {
            parent->calcVisualDensity(radiusDensityEstimation);
        }
        else {
            double volume = radius * radius * radius;
            if(volume == 0 || mass == 0) return;
            double density = mass / volume;

            if(density == 0 || density == INFINITY) return;

            size_t numChildren = childParticles.size(); // Assuming childParticles is a std::vector
            for (size_t i = 0; i < numChildren; ++i) {
                Particle* currentParticle = childParticles[i];
                
                // Check if currentParticle is valid
                if (reinterpret_cast<std::uintptr_t>(currentParticle) < 0x100000) {
                    // Handle invalid particle address
                    continue;
                }

                // Proceed with calculations using currentParticle
                childParticles[i]->visualDensity = density; // Crash occurs here
            }
        }
    }
}

