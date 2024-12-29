#ifndef SIMULATION_SETTINGS_H_
#define SIMULATION_SETTINGS_H_

#include <stdint.h>
#include <string>

// Contains all settings necessary to run the simulation.
struct SimulationSettings
{
    std::string path;
    // Formerly wrongly known as speciation_rate_per_habitat
    uint64_t speciation_rate_per_population;
    double mean_bodymass_ratio_predator_prey;
    double initial_dispersal_rate;
    double dispersal_variance;
    double zero_crossing;
    double min_feeding_range;
    double max_feeding_range;
    size_t grid_length;

    // Default constructor
    SimulationSettings(
        uint64_t speciation_rate_per_population,
        double mean_bodymass_ratio_predator_prey,
        double initial_dispersal_rate,
        double dispersal_variance,
        double zero_crossing,
        double min_feeding_range,
        double max_feeding_range,
        size_t grid_length
        ) : speciation_rate_per_population(speciation_rate_per_population),
        mean_bodymass_ratio_predator_prey(mean_bodymass_ratio_predator_prey),
        initial_dispersal_rate(initial_dispersal_rate),
        dispersal_variance(dispersal_variance),
        zero_crossing(zero_crossing),
        min_feeding_range(min_feeding_range),
        max_feeding_range(max_feeding_range),
        grid_length(grid_length)
    {

    }

    // Returns simulation setting defaults
    static const SimulationSettings DEFAULT()
    {
        return SimulationSettings(
            1000000, //speciation_rate_per_population
            2.0, //mean_bodymass_ratio_predator_prey
            10.0, //initial_dispersal_rate
            0.3, //dispersal_variance
            1e5, //zero_crossing
            0.2, //min_feeding_range
            2.0, //max_feeding_range
            5
            );
    }

    size_t maximum_species_number()
    {
        return 30 * grid_size() + 1001;
    }

    size_t grid_size()
    {
        return grid_length * grid_length;
    }
};

#endif