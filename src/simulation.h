#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <vector>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "base_settings.h"
#include "simulation_settings.h"
#include "foodweb.h"
#include "species.h"

class Simulation
{
public:
    // public constructor - will return a non-functional simulation
    Simulation();
    /// @brief formerly known as init - will return a new simulation from the given settings
    /// @param base_settings 
    /// @param settings 
    Simulation(BaseSettings base_settings, SimulationSettings settings);
    /// @brief will return a simulation loaded from the file
    /// @param base_settings 
    /// @param path path to load from
    Simulation(BaseSettings base_settings, std::string path);
    void run();

    int m_result;

private:
    // formerly known as basic init
    Simulation(SimulationSettings settings, double speciations_per_patch);

    SimulationSettings m_settings;
    /// @brief starts a speciation
    bool handle_speciation(size_t &x, size_t &y);
    /// @brief Finds a web to speciate from
    /// @param[out] target_x x coordinate of the foodweb to speciate from
    /// @param[out] target_y y coordinate of the foodweb to speciate from
    /// @return true if a web has been found, false otherwise
    bool find_web_for_speciation(size_t &target_x, size_t &target_y);
    uint64_t find_species_for_speciation(size_t target_x, size_t target_y);
    Species* speciate(Species *parent);
    /// @brief starts a dispersal
    void handle_dispersal();

    /// @brief Calculates predatorial strength
    /// @return 
    double calculate_predator_strength(double dispersal_rate);
    /// @brief returns a random value between 0.0 and 1.0 (exclusive)
    /// @return 
    double random_value();
    bool die(size_t x, size_t y);

    gsl_rng *m_generator;

    double m_t;
    uint64_t m_speciation_rate_per_population;
    /// @brief This will initially equal to settings.initial_dispersal_rate * settings.speciation_rate_per_population
    double m_initial_dispersal_rate;
    uint64_t m_total_dispersal_rate;
    double m_zero_crossing;
    // current sum of populations
    // AKA P
    uint64_t m_population_count;


    double m_speciations_per_patch;

    Species **m_species;
    /// @brief list of how many of a species exist in patches
    size_t *m_species_count;
    std::vector<size_t> m_free_indices;
    /// @brief How many different species are in our m_species array
    /// AKA S
    size_t m_number_of_living_species;

    Foodweb ***m_foodwebs;

private:
    enum ErrorCodes
    {
        WebNotFound = 1
        , FoodwebFull
        , TooManySpecies
    };
};

#endif