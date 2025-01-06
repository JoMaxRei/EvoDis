#include "simulation.h"

#include <algorithm>
#include <iostream>
// #include <sstream>
#include "easylogging++.h"
#include "exception.h"
#include "foodweb_cache.h"

Simulation::Simulation() : m_result(-1), m_settings(SimulationSettings::DEFAULT())
{
}

Simulation::Simulation(BaseSettings base_settings, SimulationSettings settings) : Simulation(settings, base_settings.speciations_per_patch)
{
    m_t = 0.0;
    LOG(INFO) << m_t << " - START";
    m_initial_dispersal_rate = settings.initial_dispersal_rate * (double)settings.speciation_rate_per_population;
    m_zero_crossing = settings.zero_crossing * (double)settings.speciation_rate_per_population;
    gsl_rng_set(m_generator, settings.seed);

    m_species_count[1] = settings.grid_size();
    m_species[1] = new Species(
        2.0,
        0.0,
        (settings.min_feeding_range + settings.max_feeding_range) / 2.0,
        calculate_predator_strength(m_initial_dispersal_rate),
        0.0,
        (uint64_t)m_initial_dispersal_rate,
        1);
    size_t next_index = m_free_indices.back();
    m_free_indices.pop_back();
    m_species[1]->set_position(next_index);
    m_total_dispersal_rate = (uint64_t)settings.grid_size() * (uint64_t)m_initial_dispersal_rate;
    for (size_t x = 0; x < settings.grid_length; x++)
    {
        for (size_t y = 0; y < settings.grid_length; y++)
        {
            m_foodwebs[x][y]->add_species(m_species[1]);
        }
    }

    m_population_count = (uint64_t)settings.grid_size();
    m_number_of_living_species += 1;

    // Calculate equilibrium for one to hash all foodwebs that are initialized the same
    // which they are
    // TODO: fix
    // FoodwebCache::calculate_equilibrium(m_foodwebs[0][0]);

    LOG(DEBUG) << m_t << " - END";
}

Simulation::Simulation(BaseSettings base_settings, std::string path) : Simulation()
{
}

void Simulation::run()
{
    while (m_result == 0 && m_t < m_speciations_per_patch)
    {
        try
        {

            double total_speciation_rate = (double)m_speciation_rate_per_population * m_population_count;
            bool event_is_speciation = (((double)m_total_dispersal_rate + total_speciation_rate) * random_value() >= (double)m_total_dispersal_rate);
            size_t x;
            size_t y;
            bool recalculate_equilibrium = false;
            if (event_is_speciation)
            {
                recalculate_equilibrium = handle_speciation(x, y);
            }
            else
            {
                handle_dispersal();
            }

            if (recalculate_equilibrium)
            {
                std::vector<size_t> dead_pops = FoodwebCache::calculate_equilibrium(m_foodwebs[x][y]);
            }

            m_t += 1.0;
        }
        catch (Exception &e)
        {
            LOG(ERROR) << e.message();
            m_result = e.error_code();
        }
    }
}

Simulation::Simulation(SimulationSettings settings, double speciations_per_patch) : m_result(0),
                                                                                    m_speciation_rate_per_population(settings.speciation_rate_per_population),
                                                                                    m_speciations_per_patch(speciations_per_patch),
                                                                                    m_settings(settings)
{
    for (size_t index = settings.maximum_species_number() - 1; index > 0; index--)
    {
        m_free_indices.push_back(index);
    }

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
    m_number_of_living_species = 1;

    // Random Number Generator initialisieren
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    // // default random number generator (so called mt19937)
    T = gsl_rng_default;
    m_generator = gsl_rng_alloc(T);
}

bool Simulation::handle_speciation(size_t &x, size_t &y)
{
    LOG(DEBUG) << m_t << " - speciating...";
    // find web
    if (!find_web_for_speciation(x, y))
    {
        throw Exception("Could not find a web for speciation", (int)ErrorCodes::WebNotFound);
    }
    if (m_foodwebs[x][y]->is_full())
    {
        // std::stringstream message;
        // message << "Foodweb x=" << x << ", y=" << y << " is full";
        throw Exception("Foodweb is full", (int)ErrorCodes::FoodwebFull);
    }
    Species *parent = m_foodwebs[x][y]->find_species_for_speciation(random_value());
    LOG(DEBUG) << "Found species " << parent->m_universal_id << " to speciate";
    Species *new_species = speciate(parent);
    if (new_species == NULL)
    {
        LOG(DEBUG) << "Speciation failed";
        return false;
    }

    if (!FoodwebCache::has_prey(m_foodwebs[x][y], new_species))
    {
        LOG(DEBUG) << "Speciation failed, species has no prey";
        delete new_species;
        return false;
    }

    if (!FoodwebCache::can_surive(m_foodwebs[x][y], new_species, m_settings))
    {
        LOG(DEBUG) << "Speciation failed, species cannot survive";
        delete new_species;
        return false;
    }

    if (m_free_indices.empty())
    {
        throw Exception("Too many species", (int)ErrorCodes::TooManySpecies);
    }
    size_t next_index = m_free_indices.back();
    m_free_indices.pop_back();

    new_species->set_position(next_index);
    m_species[next_index] = new_species;
    m_species_count[next_index] = 1;

    m_population_count += 1;
    m_number_of_living_species += 1;
    m_total_dispersal_rate += parent->m_dispersal_rate;
    return true;
}

bool Simulation::find_web_for_speciation(size_t &target_x, size_t &target_y)
{
    size_t sum1 = (size_t)((double)m_population_count * random_value());
    size_t sum2 = 0;
    LOG(DEBUG) << "sum1 is " << sum1;
    for (size_t x = 0; x < m_settings.grid_length; x++)
    {
        for (size_t y = 0; y < m_settings.grid_length; y++)
        {
            LOG(DEBUG) << "Foodweb dimension for x=" << x << ", y=" << y << " is " << m_foodwebs[x][y]->get_dimension();
            sum2 += m_foodwebs[x][y]->get_dimension() - 1;
            LOG(DEBUG) << "sum2 is " << sum2;
            if (sum1 < sum2)
            {
                target_x = x;
                target_y = y;
                return true;
            }
        }
    }
    return false;
}

Species *Simulation::speciate(Species *parent)
{
    double new_bodymass = parent->m_bodymass + 2.0 * log10(5.0) * (random_value() - 0.5);
    if (new_bodymass <= 0.0)
    {
        return NULL;
    }

    double new_feeding_center;
    double new_feeding_range;

    do
    {
        do
        {
            new_feeding_center = new_bodymass - m_settings.mean_bodymass_ratio_predator_prey + gsl_ran_gaussian(m_generator, 1);
        } while (new_feeding_center - m_settings.max_feeding_range > new_bodymass);

        new_feeding_range = m_settings.min_feeding_range + random_value() * (m_settings.max_feeding_range - m_settings.min_feeding_range);
    } while (new_feeding_center - new_feeding_range > new_bodymass);

    double new_dispersal_rate_min = log2((double)parent->m_dispersal_rate / (1.0 + m_settings.dispersal_variance));
    double new_dispersal_rate_max = log2((double)parent->m_dispersal_rate * (1.0 + m_settings.dispersal_variance));
    double new_dispersal_rate_diff = new_dispersal_rate_max - new_dispersal_rate_min;
    double new_dispersal_rate = round(pow(2, new_dispersal_rate_min + random_value() * new_dispersal_rate_diff));
    double new_predator_strength = calculate_predator_strength(new_dispersal_rate);

    return new Species(
        new_bodymass, new_feeding_center, new_feeding_range, new_predator_strength, m_t, (uint64_t)new_dispersal_rate);
}

void Simulation::handle_dispersal()
{
    LOG(DEBUG) << m_t << " - dispersing...";
}

double Simulation::calculate_predator_strength(double dispersal_rate)
{
    double value = log10(dispersal_rate / m_zero_crossing) / log10(m_initial_dispersal_rate / m_zero_crossing);
    // return max(0.0, min(1.0, value));
    return std::clamp(value, 0.0, 1.0);
}

double Simulation::random_value()
{
    return gsl_rng_uniform(m_generator);
}
