#include "foodweb.h"

Foodweb::Foodweb(Species *resource, size_t x, size_t y)
{
    m_x = x;
    m_y = y;

    m_species = new Species *[MAX_DIM];
    m_fitness = new double[MAX_DIM];
    m_species[0] = resource;
    m_species_count = 1;
}

size_t Foodweb::get_dimension() const
{
    return m_species_count;
}

int64_t Foodweb::add_species(Species *species)
{
    if (m_species_count >= MAX_DIM)
    {
        return -1;
    }

    size_t position_in_foodweb;
    for(size_t i = m_species_count - 1; i >= 0; i--)
    {
        if (species->m_bodymass < m_species[i]->m_bodymass)
        {
            m_species[i + 1] = m_species[i];
            // appearance_time[i+1] = appearance_time[i];
        }
        else
        {
            m_species[i + 1] = species;
            // appearance_time[i+1] = time;
            position_in_foodweb = i + 1;
            break;
        }
    }

    m_species_count += 1;

    m_local_dispersal_rate += species->m_dispersal_rate;
    return (int64_t)position_in_foodweb;
}