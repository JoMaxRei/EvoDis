#include "foodweb.h"


Foodweb::Foodweb(Species *resource, int64_t id)
{
    m_id = id;

    m_species = new Species*[MAX_DIM];
    m_fitness = new double[MAX_DIM];
    m_species[0] = resource;
    m_species_count = 1;
}

int64_t Foodweb::add_species(Species *species)
{
    return -1;
}