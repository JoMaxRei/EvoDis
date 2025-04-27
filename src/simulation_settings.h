#ifndef SIMULATION_SETTINGS_H_
#define SIMULATION_SETTINGS_H_

#include <cstddef>
#include <stdint.h>
#include <string>

// Contains all settings necessary to run the simulation.
struct SimulationSettings
{
    /// @brief Formerly wrongly known as speciation_rate_per_habitat
    uint64_t speciation_rate_per_population;
    /// @brief Mean difference in bodymass between predator and prey.
    double mean_bodymass_ratio_predator_prey;
    double initial_dispersal_rate;
    double dispersal_variance;
    double zero_crossing;
    double min_feeding_range;
    double max_feeding_range;
    /// @brief AKA n OR habitat grid length
    size_t grid_length;
    /// @brief AKA l
    size_t dispersal_range;
    bool periodic_boundary_conditions;
    unsigned long seed;
    /// @brief AKA p
    double predation_parameter;
    /// @brief AKA E0
    double base_gain;
    /// @brief AKA ξ
    double competition_parameter;
    /// @brief AKA x
    double trophic_loss;

    // Default constructor
    SimulationSettings(
        uint64_t speciation_rate_per_population
        , double mean_bodymass_ratio_predator_prey
        , double initial_dispersal_rate
        , double dispersal_variance
        , double zero_crossing
        , double min_feeding_range
        , double max_feeding_range
        , size_t grid_length
        , size_t dispersal_range
        , bool periodic_boundary_conditions
        , unsigned long seed
        , double predation_parameter
        , double base_gain
        , double competition_parameter
        , double trophic_loss
        ) :
        speciation_rate_per_population(speciation_rate_per_population)
        , mean_bodymass_ratio_predator_prey(mean_bodymass_ratio_predator_prey)
        , initial_dispersal_rate(initial_dispersal_rate)
        , dispersal_variance(dispersal_variance)
        , zero_crossing(zero_crossing)
        , min_feeding_range(min_feeding_range)
        , max_feeding_range(max_feeding_range)
        , grid_length(grid_length)
        , dispersal_range(dispersal_range)
        , periodic_boundary_conditions(periodic_boundary_conditions)
        , seed(seed)
        , predation_parameter(predation_parameter)
        , base_gain(base_gain)
        , competition_parameter(competition_parameter)
        , trophic_loss(trophic_loss)
    {

    }

    // Returns simulation setting defaults
    static const SimulationSettings DEFAULT()
    {
        return SimulationSettings(
            1000000     // speciation_rate_per_population
            , 2.0       // mean_bodymass_ratio_predator_prey
            , 10.0      // initial_dispersal_rate
            , 0.3       // dispersal_variance
            , 1e5       // zero_crossing
            , 0.2       // min_feeding_range
            , 2.0       // max_feeding_range
            , 5         // grid length
            , 1         // dispersal range
            , true      // periodic boundary conditions
            , 0         // seed
            , 1.0       // predation parameter
            , 1000.0    // base gain
            , 3.0       // competition parameter
            , 0.4       // trophic loss
            );
    }

    /// @brief AKA n * n = n ^ 2
    ///
    /// AKA grid size
    /// @return grid_length squared
    size_t number_of_habitats()
    {
        return grid_length * grid_length;
    }

    size_t maximum_species_number()
    {
        return 30 * number_of_habitats() + 1001;
    }
    
};

#endif