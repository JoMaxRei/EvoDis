#ifndef FOODWEB_CACHE_H_
#define FOODWEB_CACHE_H_

#include <vector>
#include "foodweb.h"
#include "species.h"
#include "simulation_settings.h"

class FoodwebCache
{
public:
    /// @brief Checks if a species has prey and can eventually survive
    /// @param target 
    /// @param species 
    /// @return true if the species can be inserted, false otherwise
    static bool has_prey(Foodweb *target, Species *new_species);

    /// @brief Checks if a species can survive in the foodweb
    ///
    /// If the species can survive in the foodweb, it is added and the foodweb has been recalculated
    /// after this method ends.
    /// Does NOT delete the new species, if it couldn't be added!
    /// @param target 
    /// @param new_species 
    /// @return true if the species can survive, false otherwise
    static bool can_surive(Foodweb *target, Species *new_species, SimulationSettings settings);

    /// @brief Calculates equilibrium of a foodweb
    ///
    /// Alters the original foodweb, if necessary
    /// @param[in,out] target 
    /// @return a vector of the global indices of species that have died
    static std::vector<size_t> calculate_equilibrium(Foodweb *target, SimulationSettings settings);
};

#endif