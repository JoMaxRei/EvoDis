#include "species.h"

Species::Species(
        double bodymass,
        double feeding_center,
        double feeding_range,
        double predator_strength,
        uint64_t dispersal_rate,
        double first_occurence,
        uint64_t universal_id
) : m_position_in_array(0),
    m_bodymass(bodymass),
    m_feeding_center(feeding_center),
    m_feeding_range(feeding_range),
    m_predator_strength(predator_strength),
    m_dispersal_rate(dispersal_rate),
    m_first_occurence(first_occurence),
    m_universal_id(universal_id)
{

}

void Species::set_position(size_t position_in_array)
{
    m_position_in_array = position_in_array;
}


void Species::update_trophic_level(double level)
{
    m_trophic_level << level;
}

double Species::get_trophic_level() const
{
    return m_trophic_level.get_mean();
}