#ifndef FOODWEB_CACHE_H_
#define FOODWEB_CACHE_H_

#include "foodweb.h"
#include "species.h"

class FoodwebCache
{
public:
    /// @brief Checks if a species can be inserted into the given foodweb
    /// @param target 
    /// @param species 
    /// @return true if the species can be inserted (or the process hasn't been hashed before), false otherwise
    static bool can_insert(Foodweb *target, Species *new_species);

    /// @brief Returns the foodweb after all calculations and eliminations have been done
    /// @param target 
    /// @return the final state of the foodweb
    static Foodweb calculate_equilibrium(Foodweb &target);
};

#endif