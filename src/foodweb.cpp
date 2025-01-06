#include "foodweb.h"

#include "easylogging++.h"

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
    for (size_t i = m_species_count - 1; i >= 0; i--)
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

Species *Foodweb::find_species_for_speciation(double random_value)
{
    double random_id = (double)(m_species_count - 1) * random_value;
    size_t local_id = (size_t)random_id + 1;
    return m_species[local_id];
}

Species *Foodweb::get_species(size_t index)
{
    if (index >= m_species_count)
    {
        LOG(ERROR) << "Trying to access element " << index << " when species array is only " << m_species_count << " long";
    }
    return m_species[index];
}

bool Foodweb::is_full() const
{
    if (m_species_count >= MAX_DIM)
    {
        return true;
    }
    return false;
}

double Foodweb::calculate_trophic_level(size_t index)
{
    return 0.0;
}

void Foodweb::calculate(SimulationSettings settings)
{
    LOG(DEBUG) << " - START";
    std::vector<std::vector<size_t>> preys (m_species_count, std::vector<size_t>{});
    std::vector<int> number_of_preys (m_species_count, 0);
    std::vector<std::vector<double>> epsilon (m_species_count, std::vector<double>{});

    calculate_feeding_relationships(preys, number_of_preys, epsilon);

    std::vector<double> sigma (m_species_count, 0.0);

    for (size_t predator_index = 1; predator_index < m_species_count; predator_index++)
    {
        for (size_t prey_index = 0; prey_index < number_of_preys[predator_index]; prey_index++)
        {
            sigma[preys[predator_index][prey_index]] += epsilon[predator_index][prey_index];
        }
    }

    std::vector<double> predation_loss (m_species_count, 0.0);
    for(size_t index = 0; index < m_species_count; index++)
    {
        predation_loss[index] = (settings.predation_parameter * sigma[index]) / (1.0 + settings.predation_parameter * sigma[index]);
    }

    std::vector<double> E (m_species_count, 0.0);
    E[0] = settings.base_gain;
    std::vector<double> E_tilde = E;

    
    for (size_t predator_index = 1; predator_index < m_species_count; predator_index++)
    {
        for (size_t prey_index = 0; prey_index < number_of_preys[predator_index]; prey_index++)
        {
            size_t current_prey_index = preys[predator_index][prey_index];
            double gain = E_tilde[current_prey_index] * predation_loss[current_prey_index] * epsilon[predator_index][prey_index] / sigma[current_prey_index];
            E_tilde[predator_index] += gain;
            E[predator_index] += gain / ( 1.0 + settings.competition_parameter * (sigma[current_prey_index] - epsilon[predator_index][prey_index]));
        }

        E_tilde[predator_index] *= settings.trophic_loss;
        E[predator_index] *= settings.trophic_loss;
    }

    for(size_t index = 1; index < m_species_count; index++)
    {
        m_fitness[index] = 1.0 / ( predation_loss[index] + 1.0 / E[index]);
        LOG(DEBUG) << "Fitness for " << index << " is " << m_fitness[index];
    }

    LOG(DEBUG) << " - END";
}

double Foodweb::get_fitness(size_t index) const
{
    return m_fitness[index];
}

void Foodweb::remove_species(size_t index)
{
    m_species_count -= 1;
    m_local_dispersal_rate -= m_species[index]->m_dispersal_rate;

    for(size_t i = index; i < m_species_count; i++)
    {
        m_species[i] = m_species[i + 1];
    }
}

void Foodweb::calculate_feeding_relationships(
    std::vector<std::vector<size_t>> &preys
    , std::vector<int> &number_of_preys
    , std::vector<std::vector<double>> &epsilon
)
{
    LOG(DEBUG) << " - START";
    for(size_t predator_index = 1; predator_index < m_species_count; predator_index++)
    {
        Species *predator = m_species[predator_index];
        for(size_t prey_index =0; prey_index < predator_index; prey_index++)
        {
            Species *prey = m_species[prey_index];

            double relative_difference = abs(prey->m_bodymass - predator->m_feeding_center) / predator->m_feeding_range;

            if (relative_difference < 1.0)
            {
                preys[predator_index].push_back(prey_index);
                number_of_preys[predator_index] += 1;

                double predation_strength = INVERTED_SQRT_HALF_PI * predator->m_predator_strength * exp(-2.0 * relative_difference * relative_difference);
                epsilon[predator_index].push_back(predation_strength);

                if (predation_strength < 0.0)
                {
                    LOG(DEBUG) << "Epsilon for " << predator_index << " is " << predation_strength;
                }
            }
        }
    }
    LOG(DEBUG) << " - END";
}

bool Foodweb::determine_dying(size_t &index)
{
    double min_fitness = 1.0;
    for(size_t i = 1; i < m_species_count; i++)
    {
        if(m_fitness[i] < min_fitness)
        {
            min_fitness = m_fitness[i];
            index = i;
        }
    }

    return min_fitness < 1.0;
}