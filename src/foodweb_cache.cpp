#include "foodweb_cache.h"

#include "easylogging++.h"

bool FoodwebCache::has_prey(Foodweb *target, Species *new_species)
{
    for (size_t index = 0; index < target->get_dimension(); index++)
    {
        Species *species = target->get_species(index);
        if (new_species->m_bodymass < species->m_bodymass)
        {
            return false;
        }
        if (abs(species->m_bodymass - new_species->m_feeding_center) < new_species->m_feeding_range)
        {
            return true;
        }
    }
    return false;
}

bool FoodwebCache::can_surive(Foodweb *target, Species *new_species, SimulationSettings settings)
{
    bool calculated = target->save_state();

    size_t new_index = target->add_species(new_species);

    // //TODO: ONLY IF INTERESTED IN SPECIES THAT DID NOT SURVIVE AT LEAST ONCE
    // double new_trophic_level = target->calculate_trophic_level(new_index);
    // new_species->update_trophic_level(new_trophic_level);

    target->calculate(settings);
    if (target->get_fitness(new_index) < 1.0)
    {
        target->remove_species(new_index);
        if (calculated)
        {
            target->restore_state();
        }

        return false;
    }

    return true;
}

std::vector<size_t> FoodwebCache::calculate_equilibrium(Foodweb *target, SimulationSettings settings)
{
    std::vector<size_t> dead_populations{};
    size_t index;
    while (true)
    {
        target->calculate(settings);
        if (!target->determine_dying(index))
        {
            break;
        }
        dead_populations.push_back(target->get_species(index)->m_position_in_array);
        target->remove_species(index);
    }

    return dead_populations;
}