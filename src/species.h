#ifndef SPECIES_H_
#define SPECIES_H_

#include <stdint.h>
#include <cstddef>
#include "eval.h"

class Species
{
public:
    Species(
        double bodymass,
        double feeding_center = 0.0,
        double feeding_range = 1.0,
        double predator_strength = 1.0,
        uint64_t dispersal_rate = 0,
        double first_occurence = 0.0,
        uint64_t universal_id = 0
        );

    double m_bodymass;
    double m_feeding_center;
    double m_feeding_range;
    /// @brief AKA u
    double m_predator_strength;
    /// @brief first time this species was present in the simulation
    double m_first_occurence;
    uint64_t m_dispersal_rate;
    /// @brief unique id of this species
    uint64_t m_universal_id;
    /// @brief index of this species in the m_species array of the simulation
    size_t m_position_in_array;

    void set_position(size_t m_position_in_array);

    /// @brief updates the trophic level of this species
    void update_trophic_level(double level);

    /// @brief returns the mean trophic level of this species over time
    double get_trophic_level() const;
    

private:
    /// @brief Sum and number of all trophic levels this species ever had after a foodweb was stable
    Eval m_trophic_level;
};

#endif