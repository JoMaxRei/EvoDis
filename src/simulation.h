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
    /// @brief performs a speciation; if succesful:
    ///
    /// - adds the species to the foodweb
    ///
    /// - recalculates that web if species has prey
    ///
    /// - updates global values (population count, total dispersal rate, number of living species, free indicies)
    /// @param[out] x x coordinate of the foodweb the speciation has happened in
    /// @param[out] y y coordinate of the foodweb the speciation has happened in
    /// @return true if the speciation was successful, false otherwise
    bool handle_speciation(size_t &x, size_t &y);
    /// @brief performs a dispersal; if succesful:
    ///
    /// - adds the species to the foodweb
    ///
    /// - recalculates that web if species has prey
    ///
    /// - updates global values (population count, total dispersal rate)
    /// @param[out] x x coordinate of the food web that was the target of the dispersal.
    /// @param[out] y y coordinate of the food web that was the target of the dispersal.
    /// @return true if the dispersal was successful, false otherwise
    bool handle_dispersal(size_t &x, size_t &y);
    /// @brief Finds a web to speciate from
    /// @param[out] target_x x coordinate of the foodweb to speciate from
    /// @param[out] target_y y coordinate of the foodweb to speciate from
    /// @return true if a web has been found, false otherwise
    bool find_web_for_speciation(size_t &target_x, size_t &target_y);

    /// @brief Finds a web to diesperse from
    /// @param[out] target_x x coordinate of the foodweb to disperse from
    /// @param[out] target_y y coordinate of the foodweb to disperse from
    /// @return true if a web has been found, false otherwise
    bool find_web_for_dispersal(size_t &target_x, size_t &target_y);
    /// @brief Finds a web to disperse to
    /// @param[out] target_x x coordinate of the foodweb to disperse to
    /// @param[out] target_y y coordinate of the foodweb to disperse to
    void find_target_web_for_dispersal(size_t &target_x, size_t &target_y);
    
    /// @brief Creates a new species by speciating from a parent species.
    /// @param parent Parent species from which the new species will be derived.
    /// @return NULL if species too small (bodymass <= 0), new species otherwise
    Species* speciate(Species *parent);


    /// @brief Calculates predatorial strength
    /// @return 
    double calculate_predator_strength(double dispersal_rate);

    /// @brief Decreases the global and specific population counter.
    ///
    /// Deletes the species if it no longer exists on any food web and decreases the species counter.
    /// @param global_index index of the species in the global species array
    void die(size_t global_index);

    /// @brief returns a random value between 0.0 and 1.0 (exclusive)
    /// @return 
    double random_value();

    gsl_rng *m_generator;

    double m_t;
    uint64_t m_speciation_rate_per_population;
    /// @brief This will initially equal to settings.initial_dispersal_rate * settings.speciation_rate_per_population
    double m_initial_dispersal_rate;
    /// @brief This will be the sum of all populations dispersal rates
    uint64_t m_total_dispersal_rate;
    double m_zero_crossing;
    /// @brief current sum of populations
    ///
    /// AKA P
    uint64_t m_population_count;


    double m_speciations_per_patch;

    Species **m_species;
    /// @brief list of how many populations of a species exist
    size_t *m_species_count;
    std::vector<size_t> m_free_indices;
    /// @brief How many different species are in our m_species array
    ///
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