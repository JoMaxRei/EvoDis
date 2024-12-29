#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "simulation_settings.h"

class Simulation
{
public:
    static Simulation create_new(SimulationSettings settings);
    static Simulation load_from_file();
    void run();

    int m_result;
};

#endif