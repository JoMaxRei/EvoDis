#include "species.h"

Species::Species(
        double bodymass,
        double feeding_center,
        double feeding_range,
        double predator_strength,
        double first_occurence,
        int64_t dispersal_rate,
        int64_t universal_id
) : m_bodymass(bodymass),
    m_feeding_center(feeding_center),
    m_feeding_range(feeding_range),
    m_predator_strength(predator_strength),
    m_first_occurence(first_occurence),
    m_dispersal_rate(dispersal_rate),
    m_universal_id(universal_id)
{

}