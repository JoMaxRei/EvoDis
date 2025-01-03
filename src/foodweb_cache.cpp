#include "foodweb_cache.h"

bool FoodwebCache::can_insert(Foodweb *target, Species *new_species)
{
    return true;
}

Foodweb FoodwebCache::calculate_equilibrium(Foodweb &target)
{
    Foodweb new_state = Foodweb(target);
    // calculate, if the state hasn't been calculated and hashed before
    // otherwise return the hashed foodweb
    return new_state;
}