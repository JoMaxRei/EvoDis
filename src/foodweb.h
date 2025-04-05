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
    // Creates a new foodweb
    // @param resource adds a non-extinguishable species to the foodweb
    // @param id the id of this new foodweb
    Foodweb(Species *resource, size_t x, size_t y);

    /// @brief returns the number of species in this foodweb, including the resource
    /// @return 
    size_t get_dimension() const;

    /// @brief Adds a species to the foodweb
    ///
    /// Also increments m_species_count and m_local_dispersal_rate
    /// @return -1 if maximum species in this web was exceeded, position of the species otherwise
    int64_t add_species(Species *species);
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
    /// @brief returns a species within this foodweb by its index
    /// @param index 
    /// @return 
    Species* get_species(size_t index);
    bool is_full() const;
    double calculate_trophic_level(size_t index);
    /// @brief Calculates the whole foodweb
    /// @param settings 
    void calculate(SimulationSettings settings);
    double get_fitness(size_t index) const;

    /// @brief Sets index to the index of the species with the lowest fitness, if below 1.0
    /// @param index 
    /// @return true, if fitness of at least one species is below 1.0, false otherwise
    bool determine_dying(size_t &index);

    // Sum of all species dispersal rates on this foodweb
    int64_t m_local_dispersal_rate;

    // x position of this foodweb
    size_t m_x;
    // y position of this foodweb
    size_t m_y;

    // maximum species that can live in a foodweb
    static const size_t MAX_DIM = 100;


private:
    /// @brief Calculates feeding relationships of this foodweb
    /// @param[in,out] preys 
    /// @param[in,out] number_of_preys 
    /// @param[in,out] epsilon 
    void calculate_feeding_relationships(
        std::vector<std::vector<size_t>> &preys
        , std::vector<int> &number_of_preys
        , std::vector<std::vector<double>> &epsilon
        );
    Species **m_species;
    // number of species on this habitat
    size_t m_species_count;
    double *m_fitness;
    
    // true if the current status of the foodweb is calculated
    bool calculated;

    static constexpr double INVERTED_SQRT_HALF_PI = 0.7978845608028654;
};

#endif