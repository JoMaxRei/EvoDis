#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "simulation_settings.h"

class Simulation
{
public:
    // public constructor - will return a non-functional simulation
    Simulation();
    static Simulation create_new(SimulationSettings settings);
    static Simulation load_from_file();
    void run();

    int m_result;

private:
    // formerly known as basic init
    Simulation(std::string path);
    gsl_rng *generator;
};

#endif