#ifndef OUTPUT_H
#define OUTPUT_H

#include <cstddef>
#include <string>
#include <vector>
#include <array>
#include <fstream>

#include "base_settings.h"
#include "simulation_settings.h"

class Output
{
public:
    /// @brief Creates a new output object
    /// @param path path to the output file
    /// @param _number_of_habitats number of habitats in the simulation (AKA n)
    /// @param initial_dispersal_rate initial dispersal rate of the simulation
    /// @param zero_crossing zero crossing value of the simulation
    Output(std::string path, size_t _number_of_habitats, double initial_dispersal_rate, double zero_crossing);

    /// @brief Access to the individual output files
    enum resfile_type {
        OUT_SETTINGS        // 0
      , OUT_HABITAT_SPECIES // 1
      , OUT_LIVING_SPECIES  // 2
      , OUT_TROPHIC_LEVELS  // 3
      , OUT_GLOBAL_INFO     // 4
      , OUT_ALIVE           // 5
      , OUT_LIFETIME        // 6
      , OUT_LTD_SLOPE       // 7
      , OUT_ABORT           // 8
      , OUT_FILE_COUNT      // 9
    };

    /// @brief Creates new output files
    /// 
    /// @return true if all files were created successfully, false otherwise
    bool create_new_files();

    /// @brief Opens all output files
    void open_files();

    /// @brief Closes all output files
    void close_files();

    /// @brief Opens a specific output file
    /// @param f file to open
    void open_file(resfile_type f);
    /// @brief Closes a specific output file
    /// @param f file to close
    void close_file(resfile_type f);

    /// @brief Suppresses the output into of a specific output file
    /// @param f file to mute
    void mute(resfile_type f);
    /// @brief Unsuppresses the output into of a specific output file
    /// @param f file to unmute
    void unmute(resfile_type f);

    /// @brief Checks if a specific output file is muted
    /// @param f file to check
    bool muted(resfile_type f);

    /// @brief Updates the lifetime bins – the general lifetime bins and the lifetime bins of the given trophic level
    /// @param lifetime Lifetime of the species
    /// @param tl_class Trophic level of the species
    void update_bins(double lifetime, size_t tl_class);

    /// @brief Prints the parameters of the simulation
    /// @param base_settings Base settings of the simulation
    /// @param settings Simulation settings of the simulation
    /// @param new_simulation true if the simulation is a new one, false if it is loaded from file
    void print_settings(resfile_type f, BaseSettings base_settings, SimulationSettings settings, bool new_simulation);

    /// @brief Prints fitness and habitat of a population
    /// @param f OUT_STEPS
    /// @param time current time of the simulation
    /// @param habitat current habitat
    /// @param universal_id species of the population
    /// @param fitness fitness of the population
    void print_line_habitat_species(resfile_type f, double time, size_t habitat, uint64_t universal_id, double fitness);

    /// @brief Prints all information about a species
    /// @param f OUT_SPECIES (FKA OUT_LIVING)
    /// @param time current time of the simulation
    /// @param universal_id universal id of the species
    /// @param first_occurence first occurence of the species
    /// @param bodymass bodymass of the species
    /// @param feeding_center feeding center of the species
    /// @param feeding_range feeding range of the species
    /// @param dispersal_rate disersal rate of the species
    /// @param predator_strength predator strength of the species
    /// @param population_count number of habitats the species populates at the current time
    /// @param mean_trophic_level mean trophic level of all populations of the species in stable foodwebs
    void print_line_living_species(resfile_type f, double time, uint64_t universal_id, double first_occurence, double bodymass, double feeding_center, double feeding_range, double dispersal_rate, double predator_strength, size_t population_count, double mean_trophic_level);

    /// @brief Prints the mean and maximum trophic level of a habitat or of all populations of the metafoodweb
    /// @param f OUT_TROPHIC_LEVELS
    /// @param time current time of the simulation
    /// @param habitat name of the habitat OR "P" for all populations of the metafoodweb
    /// @param dimension number of populations in the habitat OR number of populations in the metafoodweb
    /// @param mean_trophic_level mean trophic level of all populations in the habitat OR of all populations in the metafoodweb
    /// @param max_trophic_level maximum trophic level of all populations in the habitat OR of all populations in the metafoodweb
    void print_line_trophic_levels(resfile_type f, double time, std::string habitat, double dimension, double mean_trophic_level, double max_trophic_level);

private:

    /// @brief Creates the names of the output files
    /// @param path path to the output file
    void create_file_names(const std::string &path);

    /// @brief Rounds a double to the nearest step
    /// @param x double to round
    /// @param step step to round to
    double floor_step(double x, double step);

    /// @brief Smallest double value that can be represented reliably
    const double machine_epsilon = 1e-12;

    /// @brief One bit per channel; 
    /// mute_flags stores which channels are muted; 
    /// Bit manipulation controls switching on/off
    size_t mute_flags;

    /// @brief Number of habitats in the simulation AKA n
    size_t number_of_habitats;

    /// @brief Contains the names of the output files
    std::vector<std::string> names;

    std::array<std::ofstream, OUT_FILE_COUNT> file;

    /// @brief Counts the frequency of lifetimes; 0: all species, 1: body mass class 1, 2: body mass class 2, ...
    std::vector<std::vector<size_t>> lifetime_bins;

    /// @brief Number (length) of lifetime bins of each body mass class
    std::vector<size_t> lifetime_bins_size;

    /// @brief Number of body mass classes
    ///
    /// tl_classes + 1 = length of the lifetime_bins vector (because the first element is all species)
    size_t tl_classes;

    /// @brief Shift of the lifetime exponent to put the species in the right bin
    double smallest_lifetime_exponent;

    /// @brief The expected smallest time interval depends via the dispersal rate on a predator strength of 0.5. (by experience)
    const double expected_mean_dispersal_strength = 0.5;

    /// @brief Number of bins between 10^x and 10^(x+1)
    const size_t inverted_binsize = 10;



};

#endif // OUTPUT_H