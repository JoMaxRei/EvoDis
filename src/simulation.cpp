#include "simulation.h"

#include <iostream>
#include "easylogging++.h"
#include "exception.h"
#include "foodweb_cache.h"

Simulation::Simulation() : m_result(-1), m_settings(SimulationSettings::DEFAULT())
{
}

Simulation::Simulation(BaseSettings base_settings, SimulationSettings settings) : Simulation(settings)
{
    m_t = 0.0;
    LOG(INFO) << m_t << " - START";
    m_initial_dispersal_rate = settings.initial_dispersal_rate * (double)settings.speciation_rate_per_population;
    m_dispersal_variance = settings.dispersal_variance;
    m_zero_crossing = settings.zero_crossing * (double)settings.speciation_rate_per_population;
    m_min_feeding_range = settings.min_feeding_range;
    m_max_feeding_range = settings.max_feeding_range;
    gsl_rng_set(m_generator, 0);


    m_species_count[1] = settings.grid_size();
    m_species[1] = new Species(
        2.0,
        0.0,
        (settings.min_feeding_range + settings.max_feeding_range) / 2.0,
        calculate_predator_strength(m_initial_dispersal_rate),
        0.0,
        m_initial_dispersal_rate,
        1);
    m_species[1]->set_position(1);
    m_total_dispersal_rate = (uint64_t)settings.grid_size() * (uint64_t)m_initial_dispersal_rate;

    for(size_t x = 0; x < settings.grid_length; x++)
    {
        for(size_t y = 0; y < settings.grid_length; y++)
        {
            m_foodwebs[x][y]->add_species(m_species[1]);
        }
    }

    // Calculate equilibrium for one to hash all foodwebs that are initialized the same
    // which they are
    FoodwebCache::calculate_equilibrium(m_foodwebs[0][0]);

    LOG(DEBUG) << m_t << " - END";
}

Simulation::Simulation(BaseSettings base_settings, std::string path) : Simulation()
{
}

void Simulation::run()
{
    while (m_result == 0 && m_t < (double)m_speciation_rate_per_population)
    {
        try
        {
            double old_t = m_t;
            double total_speciation_rate = (double)m_speciation_rate_per_population * m_population_count;
            m_t += 1.0;

            bool event_is_speciation = (((double)m_total_dispersal_rate + total_speciation_rate) * gsl_rng_uniform(m_generator) >= (double)m_total_dispersal_rate);
            if (event_is_speciation)
            {
                speciate();
            }
            else
            {
                disperse();
            }
        }
        catch (Exception &e)
        {
            LOG(ERROR) << e.message();
            m_result = e.error_code();
        }
    }
}

Simulation::Simulation(SimulationSettings settings) : m_result(0), m_speciation_rate_per_population(settings.speciation_rate_per_population), m_settings(settings)
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

    // Random Number Generator initialisieren
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    // // default random number generator (so called mt19937)
    T = gsl_rng_default;
    m_generator = gsl_rng_alloc(T);
}

void Simulation::speciate()
{
    LOG(DEBUG) << m_t << " - speciating...";
    size_t x = 0;
    size_t y = 0;
    // find web
    if (!find_web_for_speciation(x, y))
    {
        throw Exception("Could not find a web for speciation", (int)ErrorCodes::WebNotFound);
    }
}

bool Simulation::find_web_for_speciation(size_t &target_x, size_t &target_y)
{
    size_t sum1 = (size_t)((double)m_population_count * gsl_rng_uniform(m_generator));
    size_t sum2 = 0;
    for (size_t x = 0; x < m_settings.grid_length; x++)
    {
        for (size_t y = 0; y < m_settings.grid_length; y++)
        {
            if (sum1 < sum2)
            {
                target_x = x;
                target_y = y;
                return true;
            }
            sum2 += m_foodwebs[x][y]->get_dimension() - 1;
        }
    }
    return false;
}

void Simulation::disperse()
{
    LOG(DEBUG) << m_t << " - dispersing...";
}

double Simulation::calculate_predator_strength(double dispersal_rate)
{
    return 0.0;
}