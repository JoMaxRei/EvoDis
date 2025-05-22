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
#include "output.h"
#include <chrono>

class Simulation
{
public:
    // public constructor - will return a non-functional simulation
    Simulation();
    /// @brief FKA init - will return a new simulation from the given settings
    /// @param base_settings 
    /// @param settings 
    Simulation(BaseSettings base_settings, SimulationSettings settings);
    /// @brief FKA load - will return a simulation loaded from the file
    /// @param base_settings 
    /// @param path path to load from
    Simulation(BaseSettings base_settings, std::string path);
    void run();

    int m_result;

    enum LogType
    {
        INTERVAL = 0,
        SAVE = 1,
        LOG_TYPE_COUNT
    };

private:
    /// @brief FKA basic init
    Simulation(SimulationSettings settings, std::string output_path, double speciations_per_patch, std::vector<std::string> muted_outputs);
    
    struct TimeComponents {
        size_t hours;
        size_t minutes;
        size_t seconds;
    };

    static TimeComponents split_time(int64_t seconds_to_split);

    void print_time(LogType type);

    inline const char* logTypeToString(LogType type);

    void create_folder(std::string path);

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

    /// @brief Finds a habitat where a speciation takes place
    /// @param[out] target_x x coordinate of the habitat to speciate from
    /// @param[out] target_y y coordinate of the habitat to speciate from
    /// @return true if a habitat has been found, false otherwise
    bool find_habitat_for_speciation(size_t &target_x, size_t &target_y);

    /// @brief Finds a habitat to diesperse from
    /// @param[out] target_x x coordinate of the habitat to disperse from
    /// @param[out] target_y y coordinate of the habitat to disperse from
    /// @return true if a habitat has been found, false otherwise
    bool find_origin_habitat_for_dispersal(size_t &target_x, size_t &target_y);

    /// @brief Finds a random neighboring habitat to disperse to
    /// @param[out] target_x x coordinate of the habitat to disperse to
    /// @param[out] target_y y coordinate of the habitat to disperse to 
    void find_target_habitat_for_dispersal(size_t &target_x, size_t &target_y);
    
    /// @brief Creates a new species by speciating from a parent species.
    /// @param parent Parent species from which the new species will be derived.
    /// @return NULL if species too small (bodymass <= 0), new species otherwise
    Species* speciate(Species *parent);

    /// @brief Calculates predatorial strength
    /// @return 
    double calculate_predator_strength(double dispersal_rate);

    /// @brief Decreases the global and specific population counter.
    ///
    /// Deletes the species and decreases the species counter if it no longer exists in any habitat.
    /// @param global_index index of the dying species in the global species array
    void die(size_t global_index);

    /// @brief returns a random value between 0.0 and 1.0 (exclusive)
    double random_value();
    /// @brief returns a random value with μ = 0 and σ = 1
    double random_normal();

    /// @brief Prints the current state of the simulation
    void print();

    /// @brief Prints fitness and habitat of all populations
    void print_steps();

    /// @brief Prints information about all currently living species
    void print_species();

    /// @brief Prints information about dispersion and dispersal of all populations and species
    void print_global_info();

    /// @brief Prints the number and the fraction of foodwebs with more species than the resource
    void print_alive_foodwebs();

    /// @brief Prints information about all trophic levels of foodwebs and of living species
    ///
    /// habitat = -1 for mean trophic levels of all populations
    ///
    /// habitat = -2 for mean trophic levels of all species
    void print_trophic_levels();

    /// @brief Prints the distribution of lifetimes of all species and for each trophic level class
    void print_lifetime_distribution();

    /// @brief Prints the slope of the lifetime distribution for all species and for each trophic level class
    void print_lifetime_distribution_slope();

    gsl_rng *m_generator;

    /// @brief Simulation time
    double m_t;

    /// @brief Counts the number of speciations that took place
    uint64_t m_speciation_counter;

    /// @brief Saving interval
    double m_save_interval;

    uint64_t m_speciation_rate_per_population;

    /// @brief This will initially equal to settings.initial_dispersal_rate * settings.speciation_rate_per_population
    double m_initial_dispersal_rate;

    /// @brief This will be the sum of all populations dispersal rates
    uint64_t m_total_dispersal_rate;

    double m_zero_crossing;

    double m_speciations_per_patch;

    Species **m_species;
    /// @brief list of how many populations of a species exist
    size_t *m_population_count;

    /// @brief current sum of all populations
    ///
    /// AKA P
    size_t m_number_of_living_populations;

    /// @brief How many different species are in our m_species array
    ///
    /// AKA S
    size_t m_number_of_living_species;

    std::vector<size_t> m_free_indices;

    Foodweb ***m_foodwebs;

    Output *m_output;

    std::chrono::steady_clock::time_point m_start_time;

    std::chrono::steady_clock::time_point m_last_log_time[LOG_TYPE_COUNT];

    int64_t m_elapsed_seconds[LOG_TYPE_COUNT];
    int64_t m_last_elapsed_seconds[LOG_TYPE_COUNT];
    double m_progress[LOG_TYPE_COUNT];
    double m_last_progress[LOG_TYPE_COUNT];

    enum ErrorCodes
    {
        WebNotFound = 1
        , FoodwebFull
        , TooManySpecies
    };
};

#endif