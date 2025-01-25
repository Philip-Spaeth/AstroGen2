#pragma once
#include <iostream>
#include "Node.h"
#include <memory>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include "Simulation.h"

class Simulation;

class Tree
{
public:
    Tree(Simulation* sim){ simulation = sim;};
    ~Tree();

    Simulation* simulation;
    Node* root;
    
    void buildTree();
    void calculateForces();
    double calcTreeWidth();
    void calcVisualDensity();
    void calcGasDensity();
};
