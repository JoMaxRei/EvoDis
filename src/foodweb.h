#ifndef FOODWEB_H_
#define FOODWEB_H_

#include <stdint.h>
#include <cstddef>
#include <vector>
#include "species.h"
#include "simulation_settings.h"
// AKA Patch
//
// Every foodweb has species  that exist in it
// Every foodweb calculates how species interact
class Foodweb
{
public:
    /// @brief Creates a new foodweb
    /// @param resource adds a non-extinguishable species to the foodweb
    Foodweb(Species *resource);

    /// @brief returns the number of species in this foodweb, including the resource
    /// @return 
    size_t get_dimension() const;

    /// @brief Adds a species to the foodweb
    ///
    /// Also increments m_species_count and m_local_dispersal_rate
    /// @return -1 if maximum species in this web was exceeded, position of the species otherwise
    size_t add_species(Species *species);
    /// @brief Removes a species from the foodweb.
    ///
    /// Also decrements m_species_count and m_local_dispersal_rate.
    /// After calling this function every entry > m_species_count is invalid.
    /// @param index 
    void remove_species(size_t index);
    /// @brief Returns a species for speciaten based on a random value
    /// @param random_value 
    /// @return 
    Species* find_species_for_speciation(double random_value);
    /// @brief Returns a species for dispersal based on a random value
    /// @param random_value 
    /// @return 
    Species* find_species_for_dispersal(double random_value);
    /// @brief Returns a species within this foodweb by its index
    /// @param index 
    /// @return 
    Species* get_species(size_t index);
    
    /// @brief Checks if the foodweb has reached its maximum size
    /// @return true if the foodweb is full, false otherwise
    bool is_full() const;

    /// @brief Calculates the whole foodweb
    /// @param settings 
    void calculate(SimulationSettings settings);

    /// @brief Calculates and updates the trophic levels of the species in this foodweb
    ///
    /// this function is only called, if the foodweb is stable
    /// @param preys Indices of the preys of each species
    /// @param number_of_preys Number of preys of each species
    void update_trophic_levels(std::vector<std::vector<size_t>> preys, std::vector<size_t> number_of_preys);
    
    // ONLY IF INTERESTED IN SPECIES THAT DID NOT SURVIVE AT LEAST ONCE
    // /// @brief Calculates the trophic level of a single species
    // /// @param index index of the species in the foodweb
    // double calculate_trophic_level(size_t index);
    
    /// @brief Returns the fitness of a population
    /// @param index position of the population in the population array of the foodweb
    /// @return fitness of the population
    double get_fitness(size_t index) const;

    /// @brief Sets index to the index of the species with the lowest fitness, if below 1.0
    /// @param index 
    /// @return true, if fitness of at least one species is below 1.0, false otherwise
    bool determine_dying(size_t &index);

    /// @brief Sum of all species dispersal rates on this foodweb
    uint64_t m_local_dispersal_rate;

    /// @brief maximum number of species that can live in a foodweb
    static const size_t MAX_DIM = 100;


private:
    /// @brief Calculates feeding relationships of this foodweb
    /// @param[in,out] preys 
    /// @param[in,out] number_of_preys 
    /// @param[in,out] epsilon 
    void calculate_feeding_relationships(
        std::vector<std::vector<size_t>> &preys
        , std::vector<size_t> &number_of_preys
        , std::vector<std::vector<double>> &epsilon
        );

    /// @brief array of species in this foodweb
    Species **m_species;

    /// @brief number of species on this habitat
    ///
    /// FKA dim
    size_t m_species_count;
    double *m_fitness;
    
    /// @brief true if the current status of the foodweb is calculated
    bool calculated;

    static constexpr double INVERTED_SQRT_HALF_PI = 0.7978845608028654;
};

#endif