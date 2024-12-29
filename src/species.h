#ifndef SPECIES_H_
#define SPECIES_H_

#include <stdint.h>

class Species
{
public:
    Species(
        double bodymass,
        double feeding_center = 0.0,
        double feeding_range = 1.0,
        double predator_strength = 1.0,
        double first_occurence = 0.0,
        int64_t dispersal_rate = 1,
        int64_t universal_id = 1
        );

    double m_bodymass;
    double m_feeding_center;
    double m_feeding_range;
    double m_predator_strength;
    double m_first_occurence;
    int64_t m_dispersal_rate;
    int64_t m_universal_id;
    int64_t m_position_in_array;

    void set_position(int64_t m_position_in_array);
};

#endif