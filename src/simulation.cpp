#include "simulation.h"

#include <iostream>

Simulation::Simulation() : m_result(-1)
{
}

Simulation Simulation::create_new(SimulationSettings settings)
{
    Simulation simulation(settings);
    simulation.m_t = 0.0;
    simulation.m_initial_dispersal_rate = (uint64_t)(settings.initial_dispersal_rate * (double)settings.speciation_rate_per_population);
    simulation.m_dispersal_variance = settings.dispersal_variance;
    simulation.m_zero_crossing = settings.zero_crossing * (double)settings.speciation_rate_per_population;
    simulation.m_min_feeding_range = settings.min_feeding_range;
    simulation.m_max_feeding_range = settings.max_feeding_range;
    return simulation;
}

Simulation Simulation::load_from_file(SimulationSettings settings)
{
    return Simulation(settings);
}

void Simulation::run()
{
    while (m_result == 0 && m_t < (double)m_speciation_rate_per_population)
    {
        double old_t = m_t;
        double total_speciation_rate = (double)m_speciation_rate_per_population * m_population_count;
        m_t += 1.0;
    }
}

Simulation::Simulation(SimulationSettings settings) : m_result(0), m_speciation_rate_per_population(settings.speciation_rate_per_population)
{
    m_species = new Species *[settings.maximum_species_number()]
    { NULL };
    m_species_count = new size_t[settings.maximum_species_number()]{0};

    m_species_count[0] = settings.grid_size();
    Species *resource = new Species(0.0, -1.0, 0.0, 0.0, 0.0, 0, 0);
    resource->set_position(0);
    m_species[0] = resource;

    m_foodwebs = new Foodweb **[settings.grid_length];
    for (size_t x = 0; x < settings.grid_length; x++)
    {
        m_foodwebs[x] = new Foodweb *[settings.grid_length];

        for (size_t y = 0; y < settings.grid_length; y++)
        {
            m_foodwebs[x][y] = new Foodweb(resource, x, y);
        }
    }
}