#ifndef SIMULATION_SETTINGS_H_
#define SIMULATION_SETTINGS_H_

#include <stdint.h>

// Contains all settings necessary to run the simulation.
struct SimulationSettings
{
    uint64_t speciation_rate_per_individuum;
    double mean_bodymass_ratio_predator_prey;

    // Default constructor
    SimulationSettings(
        uint64_t speciation_rate_per_individuum,
        double mean_bodymass_ratio_predator_prey
        ) : speciation_rate_per_individuum(speciation_rate_per_individuum),
        mean_bodymass_ratio_predator_prey(mean_bodymass_ratio_predator_prey)
    {

    }

    // Returns simulation setting defaults
    static const SimulationSettings DEFAULT()
    {
        return SimulationSettings(
            1000000,
            2.0
            );
    }
};

#endif