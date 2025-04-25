#include "foodweb.h"
#include <cmath>

#include "easylogging++.h"

Foodweb::Foodweb(Species *resource)
{
    m_species = new Species *[MAX_DIM];
    m_fitness = new double[MAX_DIM];
    m_species[0] = resource;
    m_species_count = 1;
    m_local_dispersal_rate = 0;
}

size_t Foodweb::get_dimension() const
{
    return m_species_count;
}

size_t Foodweb::add_species(Species *species)
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

    calculated = false;

    return position_in_foodweb;
}

void Foodweb::remove_species(size_t index)
{
    m_species_count -= 1;
    m_local_dispersal_rate -= m_species[index]->m_dispersal_rate;

    for (size_t i = index; i < m_species_count; i++)
    {
        m_species[i] = m_species[i + 1];
    }

    calculated = false;
}

Species *Foodweb::find_species_for_speciation(double random_value)
{
    double random_id = static_cast<double>(m_species_count - 1) * random_value;
    // LOG(DEBUG) << "random_id is " << random_id;
    size_t local_id = static_cast<size_t>(random_id) + 1;
    // LOG(DEBUG) << "Found local species " << local_id << " to speciate";
    // LOG(DEBUG) << "First occurence of this species is: " << m_species[local_id]->m_first_occurence;
    return m_species[local_id];
}

Species *Foodweb::find_species_for_dispersal(double random_value)
{
    uint64_t sum1 = static_cast<uint64_t>(static_cast<double>(m_local_dispersal_rate) * random_value);
    // LOG(DEBUG) << "sum1 is " << sum1;
    uint64_t sum2 = 0.0;
    for (size_t i = 1; i < m_species_count; i++)
    {
        // LOG(DEBUG) << "Dispersal rate for local species " << i << " is " << m_species[i]->m_dispersal_rate;
        sum2 += m_species[i]->m_dispersal_rate;
        // LOG(DEBUG) << "sum2 is " << sum2;
        if (sum1 < sum2)
        {
            // LOG(DEBUG) << "Found local species " << i << " to speciate";
            // LOG(DEBUG) << "First occurence of this species is: " << m_species[i]->m_first_occurence;
            return m_species[i];
        }
    }
    // LOG(DEBUG) << "No species found for dispersal";
    return NULL;
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
    // LOG(DEBUG) << " - START";

    if (calculated)
    {
        // LOG(DEBUG) << "already calculated - END";
        return;
    }

    std::vector<std::vector<size_t>> preys(m_species_count, std::vector<size_t>{});
    std::vector<size_t> number_of_preys(m_species_count, 0);
    std::vector<std::vector<double>> epsilon(m_species_count, std::vector<double>{});

    calculate_feeding_relationships(preys, number_of_preys, epsilon);

    std::vector<double> sigma(m_species_count, 0.0);

    for (size_t predator_index = 1; predator_index < m_species_count; predator_index++)
    {
        for (size_t prey_index = 0; prey_index < number_of_preys[predator_index]; prey_index++)
        {
            sigma[preys[predator_index][prey_index]] += epsilon[predator_index][prey_index];
        }
    }

    std::vector<double> predation_loss(m_species_count, 0.0);
    for (size_t index = 0; index < m_species_count; index++)
    {
        predation_loss[index] = (settings.predation_parameter * sigma[index]) / (1.0 + settings.predation_parameter * sigma[index]);
    }

    std::vector<double> E(m_species_count, 0.0);
    E[0] = settings.base_gain;
    std::vector<double> E_tilde = E;

    for (size_t predator_index = 1; predator_index < m_species_count; predator_index++)
    {
        for (size_t prey_index = 0; prey_index < number_of_preys[predator_index]; prey_index++)
        {
            size_t current_prey_index = preys[predator_index][prey_index];
            double gain = E_tilde[current_prey_index] * predation_loss[current_prey_index] * epsilon[predator_index][prey_index] / sigma[current_prey_index];
            E_tilde[predator_index] += gain;
            E[predator_index] += gain / (1.0 + settings.competition_parameter * (sigma[current_prey_index] - epsilon[predator_index][prey_index]));
        }

        E_tilde[predator_index] *= settings.trophic_loss;
        E[predator_index] *= settings.trophic_loss;
    }

    bool staple = true;

    for (size_t index = 1; index < m_species_count; index++)
    {
        m_fitness[index] = 1.0 / (predation_loss[index] + 1.0 / E[index]);
        // LOG(DEBUG) << "Fitness for " << index << " is " << m_fitness[index];
        if (m_fitness[index] < 1.0)
        {
            staple = false;
        }
    }

    // LOG(DEBUG) << "Staple is " << staple;

    if (staple)
    {
        update_trophic_levels(preys, number_of_preys);
    }

    calculated = true;

    // LOG(DEBUG) << " - END";
}

void Foodweb::update_trophic_levels(std::vector<std::vector<size_t>> preys, std::vector<size_t> number_of_preys)
{
    // std::vector<double> TL_shortest_Path(m_species_count, 0.0);
    
    std::vector<double> TL_mean_preys(m_species_count, 0.0);

    for (size_t j = 1; j < m_species_count; j++)
    {
        // double min = TL_shortest_Path[preys[j][0]];
        double mean = TL_mean_preys[preys[j][0]];

        // starting from 1 because 0 is already set above
        for (size_t i = 1; i < number_of_preys[j]; i++)
        {
            // if (TL_shortest_Path[preys[j][i]] < min)
            //     min = TL_shortest_Path[preys[j][i]];
            mean += TL_mean_preys[preys[j][i]];
        }

        mean /= number_of_preys[j];

        // TL_shortest_Path[j] = min + 1.0;
        TL_mean_preys[j] = mean + 1.0;

        // species[j]->update_trophic_level(TL_shortest_Path[j]);
        m_species[j]->update_trophic_level(TL_mean_preys[j]);
    }
}

double Foodweb::get_fitness(size_t index) const
{
    return m_fitness[index];
}

void Foodweb::calculate_feeding_relationships(
    std::vector<std::vector<size_t>> &preys, std::vector<size_t> &number_of_preys, std::vector<std::vector<double>> &epsilon)
{
    // LOG(DEBUG) << " - START";
    for (size_t predator_index = 1; predator_index < m_species_count; predator_index++)
    {
        Species *predator = m_species[predator_index];
        for (size_t prey_index = 0; prey_index < predator_index; prey_index++)
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
                    // LOG(DEBUG) << "Epsilon for " << predator_index << " is " << predation_strength;
                }
            }
        }
    }
    // LOG(DEBUG) << " - END";
}

bool Foodweb::determine_dying(size_t &index)
{
    double min_fitness = 1.0;
    for (size_t i = 1; i < m_species_count; i++)
    {
        if (m_fitness[i] < min_fitness)
        {
            min_fitness = m_fitness[i];
            index = i;
        }
    }

    return min_fitness < 1.0;
}