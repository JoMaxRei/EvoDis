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

Simulation::Simulation(BaseSettings base_settings, SimulationSettings settings) : Simulation(settings, base_settings.output_path, base_settings.speciations_per_patch, base_settings.muted_outputs)
{
    m_t = 0.0;
    LOG(INFO) << "New simulation";

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
        2.0,                                                             // bodymass
        0.0,                                                             // feeding center
        (settings.min_feeding_range + settings.max_feeding_range) / 2.0, // feeding range
        calculate_predator_strength(m_initial_dispersal_rate),           // predator strength
        static_cast<uint64_t>(m_initial_dispersal_rate),                 // dispersal rate
        0.0,                                                             // first occurence
        1                                                                // universal id
    );

    m_speciation_counter = 1;
    m_successful_speciation_counter = 1;
    m_failed_speciation_counter = 0;
    m_successful_disperal_counter = 0;
    m_failed_disperal_counter = 0;

    size_t next_index = m_free_indices.back();
    m_free_indices.pop_back();
    m_species[1]->set_position(next_index);
    m_total_dispersal_rate = static_cast<uint64_t>(settings.number_of_habitats()) * static_cast<uint64_t>(m_initial_dispersal_rate);
    for (size_t x = 0; x < settings.grid_length; x++)
    {
        for (size_t y = 0; y < settings.grid_length; y++)
        {
            m_foodwebs[x][y]->add_species(m_species[1], 0.0);
            m_foodwebs[x][y]->calculate(settings);
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

Simulation::Simulation(BaseSettings base_settings, std::string path) : Simulation()
{
}

Simulation::Simulation(SimulationSettings settings, std::string output_path, double speciations_per_patch, std::vector<std::string> muted_outputs) : m_result(0),
                                                                                                                                                     m_speciation_rate_per_population(settings.speciation_rate_per_population),
                                                                                                                                                     m_speciations_per_patch(speciations_per_patch),
                                                                                                                                                     m_settings(settings)
{
    LOG(INFO) << "Preparing simulation";

    for (size_t index = settings.maximum_species_number() - 1; index > 0; index--)
    {
        m_free_indices.push_back(index);
    }

    m_species = new Species *[settings.maximum_species_number()]
    { NULL };
    m_population_count = new size_t[settings.maximum_species_number()]{0};

    m_population_count[0] = settings.number_of_habitats();
    Species *resource = new Species(
        0.0,  // bodymass
        -1.0, // feeding center
        0.0,  // feeding range
        0.0,  // predator strength
        0,    // dispersal rate
        0.0,  // first occurence
        0     // universal id
    );
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

    for (const auto &name : muted_outputs)
    {
        if (name == "settings")
        {
            m_output->mute(Output::OUT_SETTINGS);
            LOG(INFO) << "Output \"settings\" is muted";
        }
        else if (name == "habitat")
        {
            m_output->mute(Output::OUT_HABITAT_SPECIES);
            LOG(INFO) << "Output \"habtitat_species\" is muted";
        }
        else if (name == "living")
        {
            m_output->mute(Output::OUT_LIVING_SPECIES);
            LOG(INFO) << "Output \"living_species\" is muted";
        }
        else if (name == "tl")
        {
            m_output->mute(Output::OUT_TROPHIC_LEVELS);
            LOG(INFO) << "Output \"trophic_levels\" is muted";
        }
        else if (name == "global")
        {
            m_output->mute(Output::OUT_GLOBAL_INFO);
            LOG(INFO) << "Output \"global_info\" is muted";
        }
        else if (name == "alive")
        {
            m_output->mute(Output::OUT_ALIVE_FOODWEBS);
            LOG(INFO) << "Output \"alive_foodwebs\" is muted";
        }
        else if (name == "ltd")
        {
            m_output->mute(Output::OUT_LIFETIME_DISTRIBUTION);
            LOG(INFO) << "Output \"lifetime_distribution\" is muted";
        }
        else if (name == "slope")
        {
            m_output->mute(Output::OUT_LTD_SLOPE);
            LOG(INFO) << "Output \"LTD_slope\" is muted";
        }
        else if (name == "abort")
        {
            m_output->mute(Output::OUT_ABORT);
            LOG(INFO) << "Output \"abort\" is muted";
        }
        else
            LOG(WARNING) << "Unknown output to mute: " << name;
    }
}

Simulation::TimeComponents Simulation::split_time(int64_t seconds_to_split)
{
    TimeComponents time;
    time.days = static_cast<size_t>(seconds_to_split / 86400);
    seconds_to_split %= 86400;

    time.hours = static_cast<size_t>(seconds_to_split / 3600);
    seconds_to_split %= 3600;

    time.minutes = static_cast<size_t>(seconds_to_split / 60);
    time.seconds = static_cast<size_t>(seconds_to_split % 60);
    return time;
}

void Simulation::print_time(LogType type)
{
    auto now = std::chrono::steady_clock::now();

    m_last_log_time[type] = now;
    m_last_progress[type] = m_progress[type];
    m_progress[type] = m_t / m_speciations_per_patch;
    m_last_elapsed_seconds[type] = m_elapsed_seconds[type];
    m_elapsed_seconds[type] = std::chrono::duration_cast<std::chrono::seconds>(now - m_start_time).count();

    double current_speed = (m_progress[type] - m_last_progress[type]) / static_cast<double>(m_elapsed_seconds[type] - m_last_elapsed_seconds[type]);

    int64_t estimated_seconds_left = static_cast<int64_t>((1 - m_progress[type]) / current_speed);

    auto time = Simulation::split_time(m_elapsed_seconds[type]);
    auto time_left = Simulation::split_time(estimated_seconds_left);
    auto time_total = Simulation::split_time(estimated_seconds_left + m_elapsed_seconds[type]);

    LOG(INFO) << logTypeToString(type) << " - [Progress]             - "
              << m_t << " / " << m_speciations_per_patch
              << " (" << 100 * m_progress[type] << "%)";
    LOG(INFO) << logTypeToString(type) << " - [Elapsed time]         - "
              << time.days << "d " << time.hours << "h " << time.minutes << "m " << time.seconds << "s";
    LOG(INFO) << logTypeToString(type) << " - [Estimated time left]  - "
              << time_left.days << "d " << time_left.hours << "h " << time_left.minutes << "m " << time_left.seconds << "s";
    LOG(INFO) << logTypeToString(type) << " - [Total estimated time] - "
              << time_total.days << "d " << time_total.hours << "h " << time_total.minutes << "m " << time_total.seconds << "s";
}

inline const char *Simulation::logTypeToString(LogType type)
{
    switch (type)
    {
    case INTERVAL:
        return "[INTERVAL]";
    case SAVE:
        return "[SAVE]    ";
    default:
        return "UNKNOWN";
    }
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

    m_start_time = std::chrono::steady_clock::now();

    std::fill(std::begin(m_last_log_time), std::end(m_last_log_time), m_start_time);

    std::fill(std::begin(m_elapsed_seconds), std::end(m_elapsed_seconds), 0);
    std::fill(std::begin(m_last_elapsed_seconds), std::end(m_last_elapsed_seconds), 0);

    std::fill(std::begin(m_progress), std::end(m_progress), 0.0);
    std::fill(std::begin(m_last_progress), std::end(m_last_progress), 0.0);

    auto log_interval = std::chrono::seconds(10);

    LOG(INFO) << "Start of simulation";

    while (m_result == 0 && m_t < m_speciations_per_patch)
    {
        try
        {

            double total_speciation_rate = static_cast<double>(m_speciation_rate_per_population) * static_cast<double>(m_number_of_living_populations);
            m_t += 1.0 / (1.0 + static_cast<double>(m_total_dispersal_rate) / total_speciation_rate) / static_cast<double>(m_settings.number_of_habitats());

            auto now = std::chrono::steady_clock::now();
            if (now - m_last_log_time[INTERVAL] >= log_interval)
            {

                print_time(INTERVAL);

                if (m_elapsed_seconds[INTERVAL] >= 120 && m_elapsed_seconds[INTERVAL] < 1800)
                    log_interval = std::chrono::seconds(60);

                if (m_elapsed_seconds[INTERVAL] >= 3600)
                    log_interval = std::chrono::seconds(1800);
            }

            if (m_t > time_of_next_save)
            {
                print_time(SAVE);
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
                if (recalculate_equilibrium)
                {
                    m_successful_speciation_counter++;
                }
                else
                {
                    m_failed_speciation_counter++;
                }
            }
            else
            {
                recalculate_equilibrium = handle_dispersal(x, y);
                if (recalculate_equilibrium)
                {
                    m_successful_disperal_counter++;
                }
                else
                {
                    m_failed_disperal_counter++;
                }
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

    if (!FoodwebCache::can_surive(m_foodwebs[x][y], new_species, m_t, m_settings)) // also adds species to foodweb
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

    for (size_t i = 1; i < m_foodwebs[x][y]->get_dimension(); i++)
    {
        if (m_foodwebs[x][y]->get_species(i)->m_universal_id == dispersing_species->m_universal_id)
        {
            // LOG(DEBUG) << "Dispersal failed, species already exists";
            // counter_inbound++;
            // if(global_species[chosen_id]->universal_id == 1)
            // {
            //     counter_inbound_S1++;
            // }
            return false;
        }

        if (m_foodwebs[x][y]->get_species(i)->m_bodymass > dispersing_species->m_bodymass)
        {
            break;
        }
    }

    if (m_foodwebs[x][y]->is_full())
    {
        throw Exception("Foodweb is full", static_cast<int>(ErrorCodes::FoodwebFull));
    }

    if (!FoodwebCache::has_prey(m_foodwebs[x][y], dispersing_species))
    {
        // LOG(DEBUG) << "Dispersal failed, species has no prey";
        return false;
    }

    if (!FoodwebCache::can_surive(m_foodwebs[x][y], dispersing_species, m_t, m_settings)) // also adds species to foodweb
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
    m_speciation_counter++;

    double new_bodymass = parent->m_bodymass + 2.0 * log10(5.0) * (random_value() - 0.5);
    if (new_bodymass <= 0.0)
    {
        return NULL;
    }

    // Option 1: Randomly choose a new feeding center and feeding range then check if feeding center plus feeding range is smaller than the bodymass and else try again
    double new_feeding_range;
    double new_feeding_center;

    do
    {
        new_feeding_range = m_settings.min_feeding_range + random_value() * (m_settings.max_feeding_range - m_settings.min_feeding_range);
        new_feeding_center = new_bodymass - m_settings.mean_bodymass_ratio_predator_prey + random_normal();
    } while (new_bodymass < new_feeding_center + new_feeding_range);

    // // Option 2: First randomly choose a new feeding range and then select a fitting feeding center
    // double new_feeding_range = m_settings.min_feeding_range + random_value() * (m_settings.max_feeding_range - m_settings.min_feeding_range);
    // double new_feeding_center;

    // do
    // {
    //     new_feeding_center = new_bodymass - m_settings.mean_bodymass_ratio_predator_prey + random_normal();
    // } while (new_feeding_center > new_bodymass - new_feeding_range);

    double new_dispersal_rate_min = log2(static_cast<double>(parent->m_dispersal_rate) / (1.0 + m_settings.dispersal_variance));
    double new_dispersal_rate_max = log2(static_cast<double>(parent->m_dispersal_rate) * (1.0 + m_settings.dispersal_variance));
    double new_dispersal_rate_diff = new_dispersal_rate_max - new_dispersal_rate_min;
    double new_dispersal_rate = round(pow(2, new_dispersal_rate_min + random_value() * new_dispersal_rate_diff));
    double new_predator_strength = calculate_predator_strength(new_dispersal_rate);

    return new Species(
        new_bodymass,                              // bodymass
        new_feeding_center,                        // feeding center
        new_feeding_range,                         // feeding range
        new_predator_strength,                     // predator strength
        static_cast<uint64_t>(new_dispersal_rate), // dispersal rate
        m_t,                                       // first occurence
        m_speciation_counter                       // universal id
    );
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
        double lifetime = m_t - m_species[global_index]->m_first_occurence;
        if (m_species[global_index]->get_trophic_level() > 0.0 || lifetime > 0.0)
        {
            m_output->update_bins(lifetime, m_species[global_index]->get_trophic_level());
        }

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

double Simulation::random_normal()
{
    return gsl_ran_gaussian(m_generator, 1);
}

void Simulation::print()
{
    LOG(INFO) << m_t << " - start printing";

    print_steps();

    print_species();

    print_trophic_levels();

    print_global_info();

    print_alive_foodwebs();

    print_lifetime_distribution();

    print_lifetime_distribution_slope();

    LOG(INFO) << m_t << " - end printing";
}

void Simulation::print_steps()
{
    if (!m_output->muted(Output::OUT_HABITAT_SPECIES))
    {
        // LOG(INFO) << "print habitat species";

        // Calculate all foodwebs. Only needed if foodweb_cache is fully implemented. Otherwise each foodweb is already calculated
        // calculate();

        m_output->open_file(Output::OUT_HABITAT_SPECIES);
        for (size_t i = 0; i < m_settings.grid_length; i++)
        {
            for (size_t j = 0; j < m_settings.grid_length; j++)
            {
                for (size_t k = 0; k < m_foodwebs[i][j]->get_dimension(); k++)
                {
                    m_output->print_line_habitat_species(Output::OUT_HABITAT_SPECIES,
                                                         m_t,
                                                         i * m_settings.grid_length + j,
                                                         m_foodwebs[i][j]->get_species(k)->m_universal_id,
                                                         m_foodwebs[i][j]->get_fitness(k));
                }
            }
        }
        // m_output->print_line(Output::OUT_STEPS, t, i*n + j, webs[i][j]->get_species(k)->universal_id);
        m_output->close_file(Output::OUT_HABITAT_SPECIES);
    }
}

void Simulation::print_species()
{
    if (!m_output->muted(Output::OUT_LIVING_SPECIES))
    {
        // LOG(INFO) << "print living species";

        size_t test_number_of_species = 0;

        m_output->open_file(Output::OUT_LIVING_SPECIES);
        for (size_t i = 0; i < m_settings.maximum_species_number(); i++)
        {
            if (m_species[i] != NULL)
            {
                test_number_of_species++;
                m_output->print_line_living_species(Output::OUT_LIVING_SPECIES,
                                                    m_t,
                                                    m_species[i]->m_universal_id,
                                                    m_species[i]->m_first_occurence,
                                                    m_species[i]->m_bodymass,
                                                    m_species[i]->m_feeding_center,
                                                    m_species[i]->m_feeding_range,
                                                    static_cast<double>(m_species[i]->m_dispersal_rate) / static_cast<double>(m_settings.speciation_rate_per_population),
                                                    m_species[i]->m_predator_strength,
                                                    m_population_count[i],
                                                    m_species[i]->get_trophic_level());
            }
        }

        if (test_number_of_species != m_number_of_living_species)
        {
            LOG(ERROR) << "Number of species in simulation: " << m_number_of_living_species << ", number of species in output: " << test_number_of_species;
        }

        m_output->close_file(Output::OUT_LIVING_SPECIES);
    }
}

void Simulation::print_trophic_levels()
{
    if (!m_output->muted(Output::OUT_TROPHIC_LEVELS))
    {
        // LOG(INFO) << "print trophic levels";

        double mean_trophic_level_populations = 0.0;
        double max_trophic_level_populations = 0.0;

        size_t test_total_dimension = 0;

        m_output->open_file(Output::OUT_TROPHIC_LEVELS);
        for (size_t i = 0; i < m_settings.grid_length; i++)
        {
            for (size_t j = 0; j < m_settings.grid_length; j++)
            {
                // It is not possible to simply average all mean trophic levels, as the different habitats contain different numbers of populations.
                mean_trophic_level_populations += static_cast<double>(m_foodwebs[i][j]->get_dimension() - 1) * m_foodwebs[i][j]->get_mean_trophic_level();
                if (max_trophic_level_populations < m_foodwebs[i][j]->get_max_trophic_level())
                {
                    max_trophic_level_populations = m_foodwebs[i][j]->get_max_trophic_level();
                }
                test_total_dimension += m_foodwebs[i][j]->get_dimension() - 1;

                m_output->print_line_trophic_levels(Output::OUT_TROPHIC_LEVELS,
                                                    m_t,
                                                    std::to_string(i * m_settings.grid_length + j),
                                                    m_foodwebs[i][j]->get_dimension(),
                                                    m_foodwebs[i][j]->get_mean_trophic_level(),
                                                    m_foodwebs[i][j]->get_max_trophic_level());
            }
        }

        m_output->print_line_trophic_levels(Output::OUT_TROPHIC_LEVELS,
                                            m_t,
                                            "P",
                                            m_number_of_living_populations,
                                            mean_trophic_level_populations / static_cast<double>(m_number_of_living_populations),
                                            max_trophic_level_populations);

        if (test_total_dimension != m_number_of_living_populations)
        {
            LOG(ERROR) << "Number of populations in simulation: " << m_number_of_living_populations << ", sum of dimensions: " << test_total_dimension;
        }

        m_output->close_file(Output::OUT_TROPHIC_LEVELS);
    }
}

void Simulation::print_global_info()
{
    if (!m_output->muted(Output::OUT_GLOBAL_INFO))
    {
        // LOG(INFO) << "print global info";

        double maximum_trophic_level = 0.0;

        for (size_t i = 1; i < m_settings.maximum_species_number(); i++)
        {
            if (m_species[i] != NULL)
            {
                if (m_species[i]->get_trophic_level() > maximum_trophic_level)
                {
                    maximum_trophic_level = m_species[i]->get_trophic_level();
                }
            }
        }

        size_t number_of_tl_bins = calc_tl_class(maximum_trophic_level);
        std::vector<double> lower_bound_of_tl_class(number_of_tl_bins + 1, 0);
        std::vector<double> upper_bound_of_tl_class(number_of_tl_bins + 1, 0);

        lower_bound_of_tl_class[0] = -1.0;
        upper_bound_of_tl_class[0] = -1.0;
        // lower_bound_of_tl_class[0] = (                                   1.0 - 0.5) / INVERTED_BINSIZE_OF_TL_CLASS;
        // upper_bound_of_tl_class[0] = (static_cast<double>(number_of_tl_bins) + 0.5) / INVERTED_BINSIZE_OF_TL_CLASS;

        for (size_t tl_class = 1; tl_class < number_of_tl_bins + 1; tl_class++)
        {
            lower_bound_of_tl_class[tl_class] = (static_cast<double>(tl_class) - 0.5) / INVERTED_BINSIZE_OF_TL_CLASS;
            upper_bound_of_tl_class[tl_class] = (static_cast<double>(tl_class) + 0.5) / INVERTED_BINSIZE_OF_TL_CLASS;
        }

        std::vector<size_t> number_of_living_populations(number_of_tl_bins + 1, 0);
        std::vector<size_t> number_of_living_species(number_of_tl_bins + 1, 0);

        std::vector<size_t> min_foodweb_size(number_of_tl_bins + 1, Foodweb::MAX_DIM);
        std::vector<double> mean_foodweb_size(number_of_tl_bins + 1, 0.0);
        std::vector<size_t> mean_foodweb_size_uint(number_of_tl_bins + 1, 0);
        std::vector<size_t> max_foodweb_size(number_of_tl_bins + 1, 0);

        std::vector<size_t> min_distribution(number_of_tl_bins + 1, m_settings.number_of_habitats());
        std::vector<double> mean_distribution(number_of_tl_bins + 1, 0.0);
        std::vector<size_t> mean_distribution_uint(number_of_tl_bins + 1, 0);
        std::vector<size_t> max_distribution(number_of_tl_bins + 1, 0);

        std::vector<double> min_dispersal_rate(number_of_tl_bins + 1, 0.0); // Wert ggf. spezifizieren
        std::vector<uint64_t> min_dispersal_rate_uint(number_of_tl_bins + 1, m_total_dispersal_rate);
        std::vector<double> mean_dispersal_rate_species(number_of_tl_bins + 1, 0.0);
        std::vector<uint64_t> mean_dispersal_rate_species_uint(number_of_tl_bins + 1, 0);
        std::vector<double> mean_dispersal_rate_populations(number_of_tl_bins + 1, 0.0);
        std::vector<uint64_t> mean_dispersal_rate_populations_uint(number_of_tl_bins + 1, 0);
        std::vector<double> max_dispersal_rate(number_of_tl_bins + 1, 0.0);
        std::vector<uint64_t> max_dispersal_rate_uint(number_of_tl_bins + 1, 0);

        std::vector<double> min_predator_strength(number_of_tl_bins + 1, 1.0);
        std::vector<double> mean_predator_strength_species(number_of_tl_bins + 1, 0.0);
        std::vector<double> mean_predator_strength_populations(number_of_tl_bins + 1, 0.0);
        std::vector<double> max_predator_strength(number_of_tl_bins + 1, 0.0);

        // size_t min_foodweb_size = Foodweb::MAX_DIM;
        // double mean_foodweb_size;
        // size_t mean_foodweb_size_uint = 0;
        // size_t max_foodweb_size = 0;
        // size_t min_distribution = m_settings.number_of_habitats();
        // double mean_distribution;
        // size_t mean_distribution_uint = 0;
        // size_t max_distribution = 0;
        // double min_dispersal_rate;
        // uint64_t min_dispersal_rate_uint = m_total_dispersal_rate;
        // double mean_dispersal_rate_species;
        // uint64_t mean_dispersal_rate_species_uint = 0;
        // double mean_dispersal_rate_populations = static_cast<double>(m_total_dispersal_rate) / static_cast<double>(m_settings.speciation_rate_per_population) / static_cast<double>(m_number_of_living_populations);
        // double max_dispersal_rate;
        // uint64_t max_dispersal_rate_uint = 0;
        // double min_predator_strength = 1.0;
        // double mean_predator_strength_species = 0.0;
        // double mean_predator_strength_populations = 0.0;
        // double max_predator_strength = 0.0;

        for (size_t i = 0; i < m_settings.grid_length; i++)
        {
            for (size_t j = 0; j < m_settings.grid_length; j++)
            {
                std::vector<size_t> number_of_populations_in_foodweb(number_of_tl_bins + 1, 0);

                for (size_t k = 1; k < m_foodwebs[i][j]->get_dimension(); k++)
                {
                    size_t tl_class = calc_tl_class(m_foodwebs[i][j]->get_species(k)->get_trophic_level());

                    number_of_living_populations[0] += 1;
                    number_of_living_populations[tl_class] += 1;

                    number_of_populations_in_foodweb[0] += 1;
                    number_of_populations_in_foodweb[tl_class] += 1;

                    mean_dispersal_rate_populations_uint[0] += m_foodwebs[i][j]->get_species(k)->m_dispersal_rate;
                    mean_dispersal_rate_populations_uint[tl_class] += m_foodwebs[i][j]->get_species(k)->m_dispersal_rate;

                    mean_predator_strength_populations[0] += m_foodwebs[i][j]->get_species(k)->m_predator_strength;
                    mean_predator_strength_populations[tl_class] += m_foodwebs[i][j]->get_species(k)->m_predator_strength;
                }

                for (size_t tl_class = 0; tl_class < number_of_tl_bins + 1; tl_class++)
                {
                    if (number_of_populations_in_foodweb[tl_class] < min_foodweb_size[tl_class])
                    {
                        min_foodweb_size[tl_class] = number_of_populations_in_foodweb[tl_class];
                    }

                    if (number_of_populations_in_foodweb[tl_class] > max_foodweb_size[tl_class])
                    {
                        max_foodweb_size[tl_class] = number_of_populations_in_foodweb[tl_class];
                    }

                    mean_foodweb_size_uint[tl_class] += number_of_populations_in_foodweb[tl_class];
                }
            }
        }

        for (size_t tl_class = 0; tl_class < number_of_tl_bins + 1; tl_class++)
        {

            mean_foodweb_size[tl_class] = static_cast<double>(mean_foodweb_size_uint[tl_class]) / static_cast<double>(m_settings.number_of_habitats());

            mean_dispersal_rate_populations[tl_class] = static_cast<double>(mean_dispersal_rate_populations_uint[tl_class]) / static_cast<double>(m_settings.speciation_rate_per_population) / static_cast<double>(number_of_living_populations[tl_class]);

            mean_predator_strength_populations[tl_class] /= static_cast<double>(number_of_living_populations[tl_class]);
        }

        for (size_t i = 1; i < m_settings.maximum_species_number(); i++)
        {
            if (m_species[i] != NULL)
            {
                size_t tl_class = calc_tl_class(m_species[i]->get_trophic_level());

                number_of_living_species[0] += 1;
                number_of_living_species[tl_class] += 1;

                mean_distribution_uint[0] += m_population_count[i];
                mean_distribution_uint[tl_class] += m_population_count[i];

                if (min_distribution[0] > m_population_count[i])
                {
                    min_distribution[0] = m_population_count[i];
                }
                if (min_distribution[tl_class] > m_population_count[i])
                {
                    min_distribution[tl_class] = m_population_count[i];
                }

                if (max_distribution[0] < m_population_count[i])
                {
                    max_distribution[0] = m_population_count[i];
                }
                if (max_distribution[tl_class] < m_population_count[i])
                {
                    max_distribution[tl_class] = m_population_count[i];
                }

                mean_dispersal_rate_species_uint[0] += m_species[i]->m_dispersal_rate;
                mean_dispersal_rate_species_uint[tl_class] += m_species[i]->m_dispersal_rate;

                if (min_dispersal_rate_uint[0] > m_species[i]->m_dispersal_rate)
                {
                    min_dispersal_rate_uint[0] = m_species[i]->m_dispersal_rate;
                }
                if (min_dispersal_rate_uint[tl_class] > m_species[i]->m_dispersal_rate)
                {
                    min_dispersal_rate_uint[tl_class] = m_species[i]->m_dispersal_rate;
                }

                if (max_dispersal_rate_uint[0] < m_species[i]->m_dispersal_rate)
                {
                    max_dispersal_rate_uint[0] = m_species[i]->m_dispersal_rate;
                }
                if (max_dispersal_rate_uint[tl_class] < m_species[i]->m_dispersal_rate)
                {
                    max_dispersal_rate_uint[tl_class] = m_species[i]->m_dispersal_rate;
                }

                mean_predator_strength_species[0] += m_species[i]->m_predator_strength;
                mean_predator_strength_species[tl_class] += m_species[i]->m_predator_strength;

                if (min_predator_strength[0] > m_species[i]->m_predator_strength)
                {
                    min_predator_strength[0] = m_species[i]->m_predator_strength;
                }
                if (min_predator_strength[tl_class] > m_species[i]->m_predator_strength)
                {
                    min_predator_strength[tl_class] = m_species[i]->m_predator_strength;
                }

                if (max_predator_strength[0] < m_species[i]->m_predator_strength)
                {
                    max_predator_strength[0] = m_species[i]->m_predator_strength;
                }
                if (max_predator_strength[tl_class] < m_species[i]->m_predator_strength)
                {
                    max_predator_strength[tl_class] = m_species[i]->m_predator_strength;
                }
            }
        }

        for (size_t tl_class = 0; tl_class < number_of_tl_bins + 1; tl_class++)
        {
            mean_distribution[tl_class] = static_cast<double>(mean_distribution_uint[tl_class]) / static_cast<double>(number_of_living_species[tl_class]);

            mean_dispersal_rate_species[tl_class] = static_cast<double>(mean_dispersal_rate_species_uint[tl_class]) / static_cast<double>(m_settings.speciation_rate_per_population) / static_cast<double>(number_of_living_species[tl_class]);

            min_dispersal_rate[tl_class] = static_cast<double>(min_dispersal_rate_uint[tl_class]) / static_cast<double>(m_settings.speciation_rate_per_population);
            max_dispersal_rate[tl_class] = static_cast<double>(max_dispersal_rate_uint[tl_class]) / static_cast<double>(m_settings.speciation_rate_per_population);

            mean_predator_strength_species[tl_class] /= static_cast<double>(number_of_living_species[tl_class]);
        }

        m_output->open_file(Output::OUT_GLOBAL_INFO);
        for (size_t tl_class = 0; tl_class < number_of_tl_bins + 1; tl_class++)
        {
            m_output->print_line_global_info(Output::OUT_GLOBAL_INFO,
                                             m_t,
                                             m_successful_speciation_counter,
                                             m_failed_speciation_counter,
                                             m_successful_disperal_counter,
                                             m_failed_disperal_counter,
                                             lower_bound_of_tl_class[tl_class],
                                             upper_bound_of_tl_class[tl_class],
                                             number_of_living_species[tl_class],
                                             number_of_living_populations[tl_class],
                                             min_foodweb_size[tl_class],
                                             mean_foodweb_size[tl_class],
                                             max_foodweb_size[tl_class],
                                             min_distribution[tl_class],
                                             mean_distribution[tl_class],
                                             max_distribution[tl_class],
                                             min_dispersal_rate[tl_class],
                                             mean_dispersal_rate_species[tl_class],
                                             mean_dispersal_rate_populations[tl_class],
                                             max_dispersal_rate[tl_class],
                                             min_predator_strength[tl_class],
                                             mean_predator_strength_species[tl_class],
                                             mean_predator_strength_populations[tl_class],
                                             max_predator_strength[tl_class]);
        }
        m_output->close_file(Output::OUT_GLOBAL_INFO);
    }
}

size_t Simulation::calc_tl_class(double trophic_level)
{
    return static_cast<size_t>(std::round(trophic_level * INVERTED_BINSIZE_OF_TL_CLASS) + MACHINE_EPSILON);
}

void Simulation::print_alive_foodwebs()
{
    if (!m_output->muted(Output::OUT_ALIVE_FOODWEBS))
    {
        // LOG(INFO) << "print alive foodwebs";

        size_t alive = 0;
        for (size_t i = 0; i < m_settings.grid_length; i++)
        {
            for (size_t j = 0; j < m_settings.grid_length; j++)
            {
                if (m_foodwebs[i][j]->get_dimension() > 1)
                {
                    alive += 1;
                }
            }
        }

        m_output->open_file(Output::OUT_ALIVE_FOODWEBS);
        m_output->print_line_alive_foodwebs(Output::OUT_ALIVE_FOODWEBS,
                                            m_t,
                                            alive,
                                            static_cast<double>(alive) / static_cast<double>(m_settings.number_of_habitats()));
        m_output->close_file(Output::OUT_ALIVE_FOODWEBS);
    }
}

void Simulation::print_lifetime_distribution()
{
    if (!m_output->muted(Output::OUT_LIFETIME_DISTRIBUTION))
    {
        // LOG(INFO) << "print lifetime distribution";

        m_output->open_file(Output::OUT_LIFETIME_DISTRIBUTION);

        m_output->print_lifetime_distribution(Output::OUT_LIFETIME_DISTRIBUTION, m_t);

        m_output->close_file(Output::OUT_LIFETIME_DISTRIBUTION);
    }
}

void Simulation::print_lifetime_distribution_slope()
{
    if (!m_output->muted(Output::OUT_LTD_SLOPE))
    {
        // LOG(INFO) << "print lifetime distribution slope";

        m_output->open_file(Output::OUT_LTD_SLOPE);

        m_output->print_LTD_slope(Output::OUT_LTD_SLOPE, m_t);

        m_output->close_file(Output::OUT_LTD_SLOPE);
    }
}