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
    enum resfile_type
    {
        OUT_SETTINGS // 0
        ,
        OUT_HABITAT_SPECIES // 1
        ,
        OUT_LIVING_SPECIES // 2
        ,
        OUT_TROPHIC_LEVELS // 3
        ,
        OUT_GLOBAL_INFO // 4
        ,
        OUT_ALIVE_FOODWEBS // 5
        ,
        OUT_LIFETIME_DISTRIBUTION // 6
        ,
        OUT_LTD_SLOPE // 7
        ,
        OUT_ABORT // 8
        ,
        // OUT_SIMULATION_TIME // 9
        // ,
        OUT_FILE_COUNT // 10
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
    /// @param trophic_level Trophic level of the species
    void update_bins(double lifetime, double ttrophic_level);

    /// @brief Prints the parameters of the simulation
    /// @param base_settings Base settings of the simulation
    /// @param settings Simulation settings of the simulation
    /// @param new_simulation true if the simulation is a new one, false if it is loaded from file
    void print_settings(resfile_type f,
                        BaseSettings base_settings,
                        SimulationSettings settings,
                        bool new_simulation);

    /// @brief Prints fitness and habitat of a population
    /// @param f OUT_STEPS
    /// @param time current time of the simulation
    /// @param habitat current habitat
    /// @param universal_id species of the population
    /// @param fitness fitness of the population
    void print_line_habitat_species(resfile_type f,
                                    double time,
                                    size_t habitat,
                                    uint64_t universal_id,
                                    double fitness);

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
    void print_line_living_species(resfile_type f,
                                   double time,
                                   uint64_t universal_id,
                                   double first_occurence,
                                   double bodymass,
                                   double feeding_center,
                                   double feeding_range,
                                   double dispersal_rate,
                                   double predator_strength,
                                   size_t population_count,
                                   double mean_trophic_level);

    /// @brief Prints mean and maximum trophic level of a habitat or of all populations of the metafoodweb
    /// @param f OUT_TROPHIC_LEVELS
    /// @param time current time of the simulation
    /// @param habitat name of the habitat OR "P" for all populations of the metafoodweb
    /// @param dimension number of populations in the habitat OR number of populations in the metafoodweb
    /// @param mean_trophic_level mean trophic level of all populations in the habitat OR of all populations in the metafoodweb
    /// @param max_trophic_level maximum trophic level of all populations in the habitat OR of all populations in the metafoodweb
    void print_line_trophic_levels(resfile_type f,
                                   double time,
                                   std::string habitat,
                                   double dimension,
                                   double mean_trophic_level,
                                   double max_trophic_level);

    /// @brief Prints global distribution and dispersal information
    /// @param f OUT_GLOBAL_INFO
    /// @param time current time of the simulation
    /// @param number_of_species number of species at the current time
    /// @param number_of_populations number of populations at the current time
    /// @param min_foodweb_size minimum number of populations over all foodwebs
    /// @param mean_foodweb_size mean number of populations over all foodwebs
    /// @param max_foodweb_size maximum number of populations over all foodwebs
    /// @param min_distribution minimum number of populations of a species
    /// @param mean_distribution mean number of populations of a species
    /// @param max_distribution maximum number of populations of a species
    /// @param min_dispersal_rate minimum dispersal rate of a species
    /// @param mean_dispersal_rate_species mean dispersal rate of a species
    /// @param mean_dispersal_rate_populations mean dispersal rate of all populations
    /// @param max_dispersal_rate maximum dispersal rate of a species
    /// @param min_predator_strength minimum predator strength of a species
    /// @param mean_predator_strength_species mean predator strength of a species
    /// @param mean_predator_strength_populations mean predator strength of all populations
    /// @param max_predator_strength maximum predator strength of a species
    void print_line_global_info(resfile_type f,
                                double time,
                                uint64_t successful_speciation_counter,
                                double percentage_of_successful_speciations,
                                uint64_t successful_disperal_counter,
                                double percentage_of_successful_dispersals,
                                double lower_bound_of_tl_class,
                                double upper_bound_of_tl_class,
                                size_t number_of_species,
                                size_t number_of_populations,
                                size_t min_foodweb_size,
                                double mean_foodweb_size,
                                size_t max_foodweb_size,
                                size_t min_distribution,
                                double mean_distribution,
                                size_t max_distribution,
                                double min_age_species,
                                double min_age_populations,
                                double mean_age_species,
                                double mean_age_populations,
                                double max_age_species,
                                double max_age_populations,
                                double min_dispersal_rate,
                                double mean_dispersal_rate_species,
                                double mean_dispersal_rate_populations,
                                double max_dispersal_rate,
                                double min_predator_strength,
                                double mean_predator_strength_species,
                                double mean_predator_strength_populations,
                                double max_predator_strength);

    /// @brief Prints the number and the fraction of foodwebs with more species than the resource
    /// @param f OUT_ALIVE_FOODWEBS
    /// @param time current time of the simulation
    /// @param number_of_alive_foodwebs number of foodwebs with more species than the resource
    /// @param fraction_of_alive_foodwebs fraction of foodwebs with more species than the resource
    void print_line_alive_foodwebs(resfile_type f,
                                   double time,
                                   size_t number_of_alive_foodwebs,
                                   double fraction_of_alive_foodwebs);

    /// @brief Prints the lifetime of all ever lived species binned and sorted by trophic level
    /// @param f OUT_LIFETIME_DISTRIBUTION
    /// @param time current time of the simulation
    void print_lifetime_distribution(resfile_type f, double time);

    /// @brief Prints the slope of the lifetime distribution of all ever lived species and also sorted by trophic level
    /// @param f OUT_LTD_SLOPE
    /// @param time current time of the simulation
    void print_LTD_slope(resfile_type f, double time);

    /// @brief Prints the current simulation time
    // void print_simulation_time(resfile_type f, double time);

private:
    /// @brief Creates the names of the output files
    /// @param path path to the output file
    void create_file_names(const std::string &path);

    /// @brief Rounds a double down to the nearest step
    /// @param x double to round
    /// @param step step to round to
    double floor_step(double x, double step);

    // /// @brief Rounds a double to the nearest step
    // /// @param x double to round
    // /// @param step step to round to
    // double round_step(double x, double step);

    /// @brief calculates the tl_class of a trophic level, depending on INVERTED_BINSIZE_OF_TL_CLASS
    /// @param trophic_level trophic level of the species
    /// @return tl_class of the species
    size_t calc_tl_class(double trophic_level);

    /// @brief Returns the trophic level interval of a tl_class
    /// @param tl_class tl_class of the interval
    /// @return [lower bound, upper bound) OR [0,∞) for tl_class 0
    std::string get_tl_interval_from_tl_class(size_t tl_class);

    /// @brief Returns the bin number the lifetime of a species falls into
    /// @param lifetime lifetime of the species
    /// @return Bin number the lifetime falls into, if lifetime is big enough (depending on smallest_lifetime_exponent), else 0
    size_t get_curr_lifetime_bin(double lifetime);

    /// @brief Updates the lifetime bin
    /// @param tl_class Trophic level class of the species OR 0 for all species
    /// @param curr_lifetime_bin Bin number the lifetime falls into
    void update_bin(size_t tl_class, size_t curr_lifetime_bin);

    /// @brief Calculates the slope of the lifetime distribution
    /// @param tl_class Trophic level class the slope is calculated for
    /// @return Slope of the lifetime distribution OR reason why it could not be calculated
    std::string calc_LTD_slope(size_t tl_class);

    /// @brief // Returns the width of a logarithmic bin depending on INVERTED_BINSIZE. Bins are defined based on powers of 10 with fractional exponent steps
    /// @param bin Bin number
    /// @return Width of the bin
    double bin_width(size_t bin);

    /// @brief Returns the midpoint position of a logarithmic bin depending on INVERTED_BINSIZE
    /// @param bin Bin number
    /// @return Midpoint position of the bin
    double bin_pos(size_t bin);

    /// @brief One bit per channel;
    ///
    /// mute_flags stores which channels are muted;
    ///
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
    static constexpr double EXPECTED_MEAN_DISPERSAL_STRENGTH = 0.5;

    /// @brief Number of bins between 10^x and 10^(x+1)
    static constexpr double INVERTED_BINSIZE = 10.0;

    /// @brief Number of tl class bins between n and n+1
    static constexpr double INVERTED_BINSIZE_OF_TL_CLASS = 1.0;

    /// @brief Smallest double value that can be represented reliably
    static constexpr double MACHINE_EPSILON = 1e-12;
};

#endif // OUTPUT_H