#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <string>
 #include <gsl/gsl_rng.h>
 #include <gsl/gsl_randist.h>
#include "simulation_settings.h"
#include "foodweb.h"
#include "species.h"

class Simulation
{
public:
    // public constructor - will return a non-functional simulation
    Simulation();
    // formerly known as init
    static Simulation create_new(SimulationSettings settings);
    static Simulation load_from_file(SimulationSettings settings);
    void run();

    int m_result;

private:
    // formerly known as basic init
    Simulation(SimulationSettings settings);

    SimulationSettings m_settings;
    void speciate();
    bool find_web_for_speciation(size_t &target_x, size_t &target_y);
    void disperse();
    gsl_rng *m_generator;

    double m_t;
    uint64_t m_speciation_rate_per_population;
    uint64_t m_initial_dispersal_rate;
    uint64_t m_total_dispersal_rate;
    double m_dispersal_variance;
    double m_zero_crossing;
    double m_min_feeding_range;
    double m_max_feeding_range;
    // current sum of populations
    // AKA P
    uint64_t m_population_count;

    Species **m_species;
    // list of how many of a species exist in patches
    size_t *m_species_count;

    Foodweb ***m_foodwebs;

private:
    enum ErrorCodes
    {
        WebNotFound = 1,
    };
};

#endif