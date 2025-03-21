#ifndef SPECIES_H_
#define SPECIES_H_

#include <stdint.h>
#include <cstddef>

class Species
{
public:
    Species(
        double bodymass,
        double feeding_center = 0.0,
        double feeding_range = 1.0,
        double predator_strength = 1.0,
        double first_occurence = 0.0,
        uint64_t dispersal_rate = 1,
        uint64_t universal_id = 0
        );

    double m_bodymass;
    double m_feeding_center;
    double m_feeding_range;
    double m_predator_strength;
    double m_first_occurence;
    uint64_t m_dispersal_rate;
    uint64_t m_universal_id;
    /// @brief index of this species in the m_species array of the simulation
    size_t m_position_in_array;

    void set_position(size_t m_position_in_array);
    void update_trophic_level(double level);
};

#endif