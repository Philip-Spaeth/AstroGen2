#include <iostream>
#include "Simulation.h"

int main()
{
    //Tittle of the Programm with information
    std::cout << std::endl << "<------------------------------------------------  Astro Genesis 2.0  ------------------------------------------------>"  << std::endl<< std::endl;

    //contrebutors and ownersion information
    std::cout << "developed by: Philip Spaeth and Kimi Sickinger  " << std::endl;
    //licence information and copy right information, licence: GNU General Public License v3.0
    std::cout << "This software is licensed under the GNU General Public License v3.0 " << std::endl<< std::endl;
    
    //The main function of the programm
    Simulation* simulation = new Simulation();
    if(simulation->init())
    {
        simulation->run();
    }
    else
    {
        std::cout << "Simulation could not be initialized " << std::endl;
    }
    delete simulation;
    return 0;
}