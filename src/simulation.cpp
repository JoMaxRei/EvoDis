#include "simulation.h"

#include <algorithm>
#include <iostream>
#include <cmath>
// #include <sstream>
#include "easylogging++.h"
#include "exception.h"
#include "foodweb_cache.h"

// defailt -> doesn't run because m_result is -1
Simulation::Simulation() : m_result(-1), m_settings(SimulationSettings::DEFAULT())
{
}

// new simulation AKA init
Simulation::Simulation(BaseSettings base_settings, SimulationSettings settings) : Simulation(settings, base_settings.output_path, base_settings.speciations_per_patch)
{
    m_t = 0.0;
    LOG(INFO) << "Start of a new simulation";

    create_folder(base_settings.output_path);
    create_folder(base_settings.output_path + "/saves");

    m_output->create_new_files();
    m_output->print_settings(Output::OUT_SETTINGS, base_settings, settings, true);

    m_save_interval = base_settings.save_interval();
    LOG(INFO) << "Simulation time: " << m_speciations_per_patch;
    LOG(INFO) << "Number of saves: " << base_settings.number_of_saves;
    LOG(INFO) << "Save interval:   " << m_save_interval;

    m_initial_dispersal_rate = settings.initial_dispersal_rate * static_cast<double>(settings.speciation_rate_per_population);
    m_zero_crossing = settings.zero_crossing * static_cast<double>(settings.speciation_rate_per_population);
    gsl_rng_set(m_generator, settings.seed);

    m_population_count[1] = settings.number_of_habitats();
    m_species[1] = new Species(
        2.0,
        0.0,
        (settings.min_feeding_range + settings.max_feeding_range) / 2.0,
        calculate_predator_strength(m_initial_dispersal_rate),
        0.0,
        static_cast<uint64_t>(m_initial_dispersal_rate),
        1);
    size_t next_index = m_free_indices.back();
    m_free_indices.pop_back();
    m_species[1]->set_position(next_index);
    m_total_dispersal_rate = static_cast<uint64_t>(settings.number_of_habitats()) * static_cast<uint64_t>(m_initial_dispersal_rate);
    for (size_t x = 0; x < settings.grid_length; x++)
    {
        for (size_t y = 0; y < settings.grid_length; y++)
        {
            m_foodwebs[x][y]->add_species(m_species[1]);
        }
    }

    m_number_of_living_populations = static_cast<size_t>(settings.number_of_habitats());
    m_number_of_living_species += 1;

    // Calculate equilibrium for one to hash all foodwebs that are initialized the same
    // which they are
    // TODO: fix
    // FoodwebCache::calculate_equilibrium(m_foodwebs[0][0]);

    // LOG(DEBUG) << m_t << " - END";
}

// AKA load
Simulation::Simulation(BaseSettings base_settings, std::string path) : Simulation()
{
}

// AKA basic_init
Simulation::Simulation(SimulationSettings settings, std::string output_path, double speciations_per_patch) : m_result(0),
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
    m_population_count = new size_t[settings.maximum_species_number()]{0};

    m_population_count[0] = settings.number_of_habitats();
    Species *resource = new Species(0.0, -1.0, 0.0, 0.0, 0.0, 0, 0);
    resource->set_position(0);
    m_species[0] = resource;

    m_foodwebs = new Foodweb **[settings.grid_length];
    for (size_t x = 0; x < settings.grid_length; x++)
    {
        m_foodwebs[x] = new Foodweb *[settings.grid_length];

        for (size_t y = 0; y < settings.grid_length; y++)
        {
            m_foodwebs[x][y] = new Foodweb(resource);
        }
    }
    m_number_of_living_species = 1;

    // initialise random number generator
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    // default random number generator (so called mt19937)
    T = gsl_rng_default;
    m_generator = gsl_rng_alloc(T);

    m_output = new Output(output_path, settings.number_of_habitats(), settings.initial_dispersal_rate, settings.zero_crossing);
}

void Simulation::create_folder(std::string path)
{
    const int dir_err = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        LOG(INFO) << "Folder already exists: " << path;
    }
}

void Simulation::run()
{
    double time_of_next_save = m_save_interval;

    while (m_result == 0 && m_t < m_speciations_per_patch)
    {
        try
        {

            double total_speciation_rate = static_cast<double>(m_speciation_rate_per_population) * static_cast<double>(m_number_of_living_populations);
            m_t += 1.0 / (1.0 + static_cast<double>(m_total_dispersal_rate) / total_speciation_rate) / static_cast<double>(m_settings.number_of_habitats());

            if (m_t > time_of_next_save)
            {
                print();
                time_of_next_save += m_save_interval;

                // if(check && (m_t > bildup_time))
                //   if(check_for_abortion(min_n_species, max_n_species, min_mean_TL, max_mean_TL, min_max_TL, max_max_TL))
                //     break;
            }

            bool event_is_speciation = ((static_cast<double>(m_total_dispersal_rate) + total_speciation_rate) * random_value() >= static_cast<double>(m_total_dispersal_rate));
            size_t x;
            size_t y;
            bool recalculate_equilibrium = false;
            if (event_is_speciation)
            {
                recalculate_equilibrium = handle_speciation(x, y);
            }
            else
            {
                recalculate_equilibrium = handle_dispersal(x, y);
            }

            if (recalculate_equilibrium)
            {
                std::vector<size_t> dead_pops = FoodwebCache::calculate_equilibrium(m_foodwebs[x][y], m_settings);
                for (size_t global_index : dead_pops)
                {
                    die(global_index);
                }
            }
        }
        catch (Exception &e)
        {
            LOG(ERROR) << e.message();
            m_result = e.error_code();
        }
    }
}

bool Simulation::handle_speciation(size_t &x, size_t &y)
{
    // LOG(DEBUG) << m_t << " - speciating...";
    // find web
    if (!find_habitat_for_speciation(x, y))
    {
        throw Exception("Could not find a web for speciation", static_cast<int>(ErrorCodes::WebNotFound));
    }
    if (m_foodwebs[x][y]->is_full())
    {
        throw Exception("Foodweb is full", static_cast<int>(ErrorCodes::FoodwebFull));
    }
    Species *parent = m_foodwebs[x][y]->find_species_for_speciation(random_value());
    // LOG(DEBUG) << "Found species " << parent->m_universal_id << " to speciate";
    // LOG(DEBUG) << "First occurence of this species is: " << parent->m_first_occurence;
    // LOG(DEBUG) << "Bodymass of this species is: " << parent->m_bodymass;
    // LOG(DEBUG) << "Feeding center of this species is: " << parent->m_feeding_center;
    // LOG(DEBUG) << "Feeding range of this species is: " << parent->m_feeding_range;
    Species *new_species = speciate(parent);
    if (new_species == NULL)
    {
        // LOG(DEBUG) << "Speciation failed";
        return false;
    }

    if (!FoodwebCache::has_prey(m_foodwebs[x][y], new_species))
    {
        // LOG(DEBUG) << "Speciation failed, species has no prey";
        delete new_species;
        return false;
    }

    if (!FoodwebCache::can_surive(m_foodwebs[x][y], new_species, m_settings)) // also adds species to foodweb
    {
        // LOG(DEBUG) << "Speciation failed, species cannot survive";
        delete new_species;
        return false;
    }

    if (m_free_indices.empty())
    {
        throw Exception("Too many species", static_cast<int>(ErrorCodes::TooManySpecies));
    }
    size_t next_index = m_free_indices.back();
    m_free_indices.pop_back();

    new_species->set_position(next_index);
    m_species[next_index] = new_species;
    m_population_count[next_index] = 1;

    m_number_of_living_populations += 1;
    m_number_of_living_species += 1;
    m_total_dispersal_rate += new_species->m_dispersal_rate;
    // LOG(DEBUG) << "Speciation succcesful";
    return true;
}

bool Simulation::handle_dispersal(size_t &x, size_t &y)
{
    // LOG(DEBUG) << m_t << " - dispersing...";
    // find web
    if (!find_origin_habitat_for_dispersal(x, y))
    {
        throw Exception("Could not find a web for dispersal", static_cast<int>(ErrorCodes::WebNotFound));
    }
    Species *dispersing_species = m_foodwebs[x][y]->find_species_for_dispersal(random_value());
    // LOG(DEBUG) << "Found species " << dispersing_species->m_universal_id << " to speciate";
    // LOG(DEBUG) << "First occurence of this species is: " << dispersing_species->m_first_occurence;
    // LOG(DEBUG) << "Bodymass of this species is: " << dispersing_species->m_bodymass;
    // LOG(DEBUG) << "Feeding center of this species is: " << dispersing_species->m_feeding_center;
    // LOG(DEBUG) << "Feeding range of this species is: " << dispersing_species->m_feeding_range;

    // Check if species exists on every foodweb
    if (m_population_count[dispersing_species->m_position_in_array] == m_settings.number_of_habitats())
    {
        return false;
    }

    find_target_habitat_for_dispersal(x, y);
    if (m_foodwebs[x][y]->is_full())
    {
        throw Exception("Foodweb is full", static_cast<int>(ErrorCodes::FoodwebFull));
    }

    if (!FoodwebCache::has_prey(m_foodwebs[x][y], dispersing_species))
    {
        // LOG(DEBUG) << "Dispersal failed, species has no prey";
        return false;
    }

    if (!FoodwebCache::can_surive(m_foodwebs[x][y], dispersing_species, m_settings)) // also adds species to foodweb
    {
        // LOG(DEBUG) << "Dispersal failed, species cannot survive";
        return false;
    }

    m_population_count[dispersing_species->m_position_in_array] += 1;

    m_number_of_living_populations += 1;
    m_total_dispersal_rate += dispersing_species->m_dispersal_rate;
    // LOG(DEBUG) << "Dispersal succcesful";
    return true;
}

bool Simulation::find_habitat_for_speciation(size_t &target_x, size_t &target_y)
{
    size_t sum1 = static_cast<size_t>(static_cast<double>(m_number_of_living_populations) * random_value());
    size_t sum2 = 0;
    // LOG(DEBUG) << "sum1 is " << sum1;
    for (size_t x = 0; x < m_settings.grid_length; x++)
    {
        for (size_t y = 0; y < m_settings.grid_length; y++)
        {
            // LOG(DEBUG) << "Foodweb dimension for x=" << x << ", y=" << y << " is " << m_foodwebs[x][y]->get_dimension();
            sum2 += m_foodwebs[x][y]->get_dimension() - 1;
            // LOG(DEBUG) << "sum2 is " << sum2;
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

bool Simulation::find_origin_habitat_for_dispersal(size_t &target_x, size_t &target_y)
{
    uint64_t sum1 = static_cast<uint64_t>(static_cast<double>(m_total_dispersal_rate) * random_value());
    uint64_t sum2 = 0;
    // LOG(DEBUG) << "sum1 is " << sum1;
    for (size_t x = 0; x < m_settings.grid_length; x++)
    {
        for (size_t y = 0; y < m_settings.grid_length; y++)
        {
            // LOG(DEBUG) << "Local dispersal rate for x=" << x << ", y=" << y << " is " << m_foodwebs[x][y]->m_local_dispersal_rate;
            sum2 += m_foodwebs[x][y]->m_local_dispersal_rate;
            // LOG(DEBUG) << "sum2 is " << sum2;
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

void Simulation::find_target_habitat_for_dispersal(size_t &x, size_t &y)
{
    // LOG(DEBUG) << "find_target_web_for_dispersal";

    size_t kernel_length = 2 * m_settings.dispersal_range + 1;
    size_t kernel_size = kernel_length * kernel_length;
    size_t source = 0;

    // case 1: periodic boundary conditions
    if (m_settings.periodic_boundary_conditions)
    {
        source = ((kernel_size + 1) / 2 + static_cast<size_t>((kernel_size - 1) * random_value())) % (kernel_size);

        x = (m_settings.dispersal_range * m_settings.grid_length + x + (source / kernel_length) - m_settings.dispersal_range) % m_settings.grid_length;
        y = (m_settings.dispersal_range * m_settings.grid_length + y + (source % kernel_length) - m_settings.dispersal_range) % m_settings.grid_length;
        return;
    }

    // case 2a: open boundary conditions, habitat at the edge
    if (x < m_settings.dispersal_range || y < m_settings.dispersal_range || x > m_settings.grid_length - m_settings.dispersal_range - 1 || y > m_settings.grid_length - m_settings.dispersal_range - 1)
    {

        size_t xmin = std::min(x, m_settings.dispersal_range);
        size_t ymin = std::min(y, m_settings.dispersal_range);
        size_t xr = xmin + std::min(m_settings.grid_length - x - 1, m_settings.dispersal_range) + 1;
        size_t yr = ymin + std::min(m_settings.grid_length - y - 1, m_settings.dispersal_range) + 1;
        source = (yr * xmin + ymin + 1 + static_cast<size_t>((xr * yr - 1) * random_value())) % (xr * yr);

        // cout << "(" << x + (source / yr) - xmin << "," << y + (source % yr) - ymin << ") -> (" << x << "," << y << ")" << endl;
        x += (source / yr) - xmin;
        y += (source % yr) - ymin;
        return;
    }

    // case 2b: open boundary conditions, habitat not at the edge
    source = ((kernel_size + 1) / 2 + static_cast<size_t>((kernel_size - 1) * random_value())) % (kernel_size);

    x += (source / kernel_length) - m_settings.dispersal_range;
    y += (source % kernel_length) - m_settings.dispersal_range;
    // Modulo isn't needed here, because the habitat is not at the edge and only valid values are generated
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

    double new_dispersal_rate_min = log2(static_cast<double>(parent->m_dispersal_rate) / (1.0 + m_settings.dispersal_variance));
    double new_dispersal_rate_max = log2(static_cast<double>(parent->m_dispersal_rate) * (1.0 + m_settings.dispersal_variance));
    double new_dispersal_rate_diff = new_dispersal_rate_max - new_dispersal_rate_min;
    double new_dispersal_rate = round(pow(2, new_dispersal_rate_min + random_value() * new_dispersal_rate_diff));
    double new_predator_strength = calculate_predator_strength(new_dispersal_rate);

    return new Species(
        new_bodymass, new_feeding_center, new_feeding_range, new_predator_strength, m_t, static_cast<uint64_t>(new_dispersal_rate));
}

void Simulation::die(size_t global_index)
{
    // LOG(DEBUG) << m_t << " - dying...";
    // LOG(DEBUG) << m_t << " - Species " << global_index << " has died";
    m_number_of_living_populations -= 1;
    m_population_count[global_index] -= 1;
    m_total_dispersal_rate -= m_species[global_index]->m_dispersal_rate;

    if (m_population_count[global_index] == 0)
    {
        // LOG(DEBUG) << "Last population of species " << global_index << " has died. Spiecies will be deleted";
        m_free_indices.push_back(global_index);
        m_species[global_index] = NULL;
        m_number_of_living_species -= 1;
    }
}

double Simulation::calculate_predator_strength(double dispersal_rate)
{
    // If this function is changed its inverse in output has to be changed too!
    double value = log10(dispersal_rate / m_zero_crossing) / log10(m_initial_dispersal_rate / m_zero_crossing);
    // return std::max(0.0, std::min(1.0, value));
    return std::clamp(value, 0.0, 1.0);
}

double Simulation::random_value()
{
    return gsl_rng_uniform(m_generator);
}

void Simulation::print()
{
    LOG(INFO) << m_t << " - print";

    print_steps();
}

void Simulation::print_steps()
{
    if (!m_output->muted(Output::OUT_STEPS))
    {
        LOG(INFO) << m_t << " - print_steps";

        // Calculate all foodwebs. Only needed if foodweb_cache is fully implemented. Otherwise each foodweb is already calculated
        // calculate();

        m_output->open_file(Output::OUT_STEPS);
        for (size_t i = 0; i < m_settings.grid_length; i++)
        {
            for (size_t j = 0; j < m_settings.grid_length; j++)
            {
                for (size_t k = 0; k < m_foodwebs[i][j]->get_dimension(); k++)
                {
                    m_output->print_line_steps(Output::OUT_STEPS, m_t, i * m_settings.grid_length + j, m_foodwebs[i][j]->get_species(k)->m_universal_id, m_foodwebs[i][j]->get_fitness(k));
                }
            }
        }
        // m_output->print_line(Output::OUT_STEPS, t, i*n + j, webs[i][j]->get_species(k)->universal_id);
        m_output->close_file(Output::OUT_STEPS);
    }
}

void Simulation::print_species()
{
    if (!m_output->muted(Output::OUT_SPECIES))
    {
        LOG(INFO) << m_t << " - print_species";

        size_t test_number_of_species = 0;

        m_output->open_file(Output::OUT_SPECIES);
        for (size_t i = 0; i < m_settings.maximum_species_number(); i++)
        {
            if (m_species[i] != NULL)
            {
                test_number_of_species++;
                m_output->print_line_species(Output::OUT_SPECIES, m_t, m_species[i]->m_universal_id, m_species[i]->m_first_occurence, m_species[i]->m_bodymass, m_species[i]->m_feeding_center, m_species[i]->m_feeding_range, m_species[i]->m_dispersal_rate, m_species[i]->m_predator_strength, m_population_count[i], m_species[i]->get_trophic_level());   
            }

        }

        if (test_number_of_species != m_number_of_living_species)
        {
            LOG(ERROR) << "Number of species in simulation: " << m_number_of_living_species << ", number of species in output: " << test_number_of_species;
        }
            
        m_output->close_file(Output::OUT_SPECIES);
    }
}