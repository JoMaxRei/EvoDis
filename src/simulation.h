#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "simulation_settings.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Simulation
{
public:
    static Simulation create_new(SimulationSettings settings);
    static Simulation load_from_file();
    void run();

    int m_result;

private:
    // formerly known as basic init
    Simulation();
    gsl_rng *generator;
};

#endif