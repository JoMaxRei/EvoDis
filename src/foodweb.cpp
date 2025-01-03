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
    return -1;
}