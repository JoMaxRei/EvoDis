#include "output.h"
#include "easylogging++.h"
#include <cmath>
#include <iomanip>

Output::Output(std::string path, size_t _number_of_habitats, double initial_dispersal_rate, double zero_crossing)
{
    create_file_names(path);

    number_of_habitats = _number_of_habitats;

    tl_classes = 0;

    lifetime_bins_size.push_back(1);

    std::vector<size_t> vec;
    lifetime_bins.push_back(vec);

    // The aim is to obtain a function that determines the associated smallest bin depending on the expected smallest time interval.
    // The smallest time interval depends on the number of habitats and the mean dispersal rate of the populations.
    // From experience, the mean dispersal rate corresponds approximately to that of a species with a predator strength of 0.5.
    double expected_mean_dispersal_rate = zero_crossing * std::pow(10, expected_mean_dispersal_strength * std::log10(initial_dispersal_rate / zero_crossing));
    double t_min = 1.0 / (1.0 + expected_mean_dispersal_rate) / static_cast<double>(number_of_habitats);

    smallest_lifetime_exponent = -floor_step(std::log10(t_min), 1.0 / inverted_binsize);
}

void Output::create_file_names(const std::string &path)
{
    static const std::array<std::string, OUT_FILE_COUNT> suffixes = {
        "settings" // 0
        ,
        "habtitat_species" // 1
        ,
        "living_species" // 2
        ,
        "trophic_levels" // 3
        ,
        "global_info" // 4
        ,
        "alive_foodwebs" // 5
        ,
        "lifetime_distribution" // 6
        ,
        "LTD_slope" // 7
        ,
        "abort" // 8
    };

    for (const auto &suffix : suffixes)
    {
        names.push_back(path + "/mc_" + suffix + ".out");
    }
}

bool Output::create_new_files()
{

    for (size_t i = 0; i < OUT_FILE_COUNT; i++)
    {
        file[i].open(names[i]);
        if (!file[i].is_open())
        {
            LOG(ERROR) << "Could not open file: " << names[i];
            return false;
        }
        switch (i)
        {
        case OUT_HABITAT_SPECIES:
            file[i] << "time" << "\t" << "habitat" << "\t" << "universal_id" << "\t" << "fitness" << std::endl;
            break;

        case OUT_LIVING_SPECIES:
            file[i] << "time" << "\t" << "universal_id" << "\t" << "first_occurence" << "\t" << "bodymass" << "\t" << "feeding_center" << "\t" << "feeding_range" << "\t" << "dispersal_rate" << "\t" << "predator_strength" << "\t" << "population_count" << "\t" << "mean_trophic_level" << std::endl;
            break;

        case OUT_TROPHIC_LEVELS:
            file[i] << "time" << "\t" << "habitat" << "\t" << "dimension" << "\t" << "mean_trophic_level" << "\t" << "max_trophic_level" << std::endl;
            break;

        case OUT_GLOBAL_INFO:
            file[i] << "time" << "\t" << "number_of_species" << "\t" << "number_of_populations" << "\t" << "min_distribution" << "\t" << "mean_distribution" << "\t" << "max_distribution" << "\t" << "min_dispersal_rate" << "\t" << "mean_dispersal_rate_species" << "\t" << "mean_dispersal_rate_populations" << "\t" << "max_dispersal_rate" << "\t" << "min_predator_strength" << "\t" << "mean_predator_strength_species" << "\t" << "mean_predator_strength_populations" << "\t" << "max_predator_strength" << std::endl;
            break;

        case OUT_ALIVE_FOODWEBS:
            file[i] << "time" << "\t" << "number_of_alive_foodwebs" << "\t" << "fraction_of_alive_foodwebs" << std::endl;
            break;

        case OUT_LIFETIME_DISTRIBUTION:
            file[i] << "time" << "\t" << "interval" << "\t" << "tl_class" << std::endl;
            file[i] << "bin content" << std::endl;
            file[i] << std::endl;
            break;

        case OUT_LTD_SLOPE:
            file[i] << "time" << "\t" << "TL Interval" << "\t" << "LTD slope" << std::endl;
            break;

        default:
            break;
        }

        file[i].close();
    }

    return true;
}

void Output::open_files()
{
    for (int i = 0; i < OUT_FILE_COUNT; i++)
    {

        // Check whether the output has been suppressed:
        if (muted(resfile_type(i)) || file[i].is_open())
            continue;
        file[i].open(names[i].c_str(), std::ios::out | std::ios::app);
    }
}

void Output::close_files()
{
    for (int i = 0; i < OUT_FILE_COUNT; i++)
    {
        // Check whether the output has been suppressed:
        if (muted(resfile_type(i)) || !file[i].is_open())
            continue;
        file[i].close();
    }
}

void Output::open_file(resfile_type f)
{
    if (muted(f) || file[f].is_open())
        return;

    file[f].open(names[f].c_str(), std::ios::out | std::ios::app);
}

void Output::close_file(resfile_type f)
{
    if (muted(f) || !file[f].is_open())
        return;

    file[f].close();
}

void Output::mute(resfile_type f)
{
    mute_flags |= (size_t(1) << static_cast<size_t>(f));
}

void Output::unmute(resfile_type f)
{
    mute_flags &= ~(size_t(1) << static_cast<size_t>(f));
}

bool Output::muted(resfile_type f)
{
    // Check whether the output has been suppressed:
    return mute_flags & (size_t(1) << static_cast<size_t>(f));
}

double Output::floor_step(double x, double step)
{
    return std::floor(x / step) * step;
}

// double Output::round_step(double x, double step)
// {
//     return std::round(x / step) * step;
// }

void Output::update_bins(double lifetime, double trophic_level)
{
    // Not necessary to update the bins if the output is muted
    if (muted(OUT_LIFETIME_DISTRIBUTION) && muted(OUT_LTD_SLOPE))
        return;

    size_t tl_class = calc_tl_class(trophic_level);

    // New species has bigger trophic level than any before
    while (tl_class > tl_classes)
    {

        tl_classes++;

        lifetime_bins_size.push_back(0);

        std::vector<size_t> vec;
        lifetime_bins.push_back(vec);
    }

    size_t curr_lifetime_bin = get_curr_lifetime_bin(lifetime);

    update_bin(0, curr_lifetime_bin);


    if (tl_class < 1) // Ressource
        return;

    update_bin(tl_class, curr_lifetime_bin);
}

size_t Output::calc_tl_class(double trophic_level)
{
    return static_cast<size_t>(std::round(trophic_level * inverted_binsize_of_tl_class) + machine_epsilon);
}

std::string Output::get_tl_interval_from_tl_class(size_t tl_class)
{
    if (tl_class == 0)
        return "[  0 ,  ∞ )";

    double lower = (static_cast<double>(tl_class) - 0.5) / inverted_binsize_of_tl_class;
    double upper = (static_cast<double>(tl_class) + 0.5) / inverted_binsize_of_tl_class;

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "[" << lower << "," << upper << ")";

    return oss.str();
}

size_t Output::get_curr_lifetime_bin(double lifetime)
{
    double curr_lifetime_bin_d = 0.0;

    // check if lifetime is 0
    if (lifetime > 0.0)
    {
        // check if lifetime is smaller than the smallest lifetime bin
        curr_lifetime_bin_d = (floor_step(std::log10(lifetime), 1. / static_cast<size_t>(inverted_binsize + machine_epsilon)) + smallest_lifetime_exponent) * static_cast<size_t>(inverted_binsize + machine_epsilon);
        if (curr_lifetime_bin_d < 0.0)
            return 0;
        else
            return static_cast<size_t>(curr_lifetime_bin_d + machine_epsilon);
    }
    else
        return 0;
}

void Output::update_bin(size_t tl_class, size_t curr_lifetime_bin)
{
    if (curr_lifetime_bin >= lifetime_bins_size[tl_class])
    {
        lifetime_bins_size[tl_class] = curr_lifetime_bin + 1;
        lifetime_bins[tl_class].resize(lifetime_bins_size[tl_class], 0);
    }
    lifetime_bins[tl_class][curr_lifetime_bin]++;
}

std::string Output::calc_LTD_slope(size_t tl_class)
{

    if (tl_class > tl_classes)
    {
        LOG(ERROR) << "tl_class: " << tl_class << " tl_classes: " << tl_classes;
        return "To high tl_class";
    }

    if (lifetime_bins_size[tl_class] < 0.5)
    {
        return "no data for this trophic level interval";
    }

    // double total = 0.0;                                          // not necessary for the logarithmic gradient, as only the x-axis is shifted
    // for(auto& a : lifetime_bins[tl_class]){total += a;}

    size_t n = 0;
    std::vector<double> x;
    std::vector<double> y;

    for (size_t i = 0; i < lifetime_bins_size[tl_class]; i++)
    {
        if (lifetime_bins[tl_class][i] > 0)
        {
            double width = bin_width(i);

            if (width > machine_epsilon)
            {
                x.push_back(std::log10(bin_pos(i)));
                y.push_back(std::log10(static_cast<double>(lifetime_bins[tl_class][i]) / width /* / total */));
                n++;
            }
        }
    }

    if (n == 0)
    {
        return "bin width was always 0 for this trophic level interval";
    }

    double avgX = 0.0;
    double avgY = 0.0;

    for (auto &a : x)
    {
        avgX += a;
    }
    for (auto &a : y)
    {
        avgY += a;
    }

    avgX /= n;
    avgY /= n;

    double numerator = 0.0;
    double denominator = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        numerator += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    if (denominator < machine_epsilon && denominator > -machine_epsilon)
    {
        return "denominator is 0 for this trophic level interval";
    }

    return std::to_string(numerator / denominator);
}

double Output::bin_width(size_t bin)
{
    // Compute the upper and lower bounds of the bin using pow(10, ...)
    // floor and ceil are used to ensure numerical stability and integer alignment.
    double width = 1. + floor(pow(10., static_cast<double>(bin) / inverted_binsize)) - ceil(pow(10., (static_cast<double>(bin) - 1.) / inverted_binsize));

    // Correct for potential off-by-one rounding at decade boundaries.
    // This occurs when the bin index is an exact multiple of inverted_binsize,
    // e.g., bin = 10 for inverted_binsize = 10 (i.e., exact powers of 10).
    if (std::abs(std::fmod(static_cast<double>(bin), inverted_binsize) - 1) < machine_epsilon)
        width--;

    return width;
}

double Output::bin_pos(size_t bin)
{
    // Estimate lower and upper bounds, then average them.
    // floor and ceil stabilize the values for nice, rounded midpoints.
    double pos = floor(pow(10., static_cast<double>(bin) / inverted_binsize)) + ceil(pow(10., (static_cast<double>(bin) - 1.) / inverted_binsize));

    // Similar correction as in bin_width: adjust at decade boundaries.
    if (std::abs(std::fmod(static_cast<double>(bin), inverted_binsize) - 1) < machine_epsilon)
        pos++;

    pos /= 2.0;

    return pos;
}

void Output::print_settings(resfile_type f,
                            BaseSettings base_settings,
                            SimulationSettings settings,
                            bool new_simulation)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    file[f] << "===============================================" << std::endl;
    if (new_simulation)
        file[f] << "New Simulation" << std::endl;
    else
        file[f] << "Loaded Simulation" << std::endl;

    file[f] << "===============================================" << std::endl;
    file[f] << std::endl;
    // if(load)
    //   file[f] << "old_path: " << old_path << std::endl;
    // file[f] << "path: " << path;
    // file[f] << std::endl;

    file[f] << "seed:                      " << settings.seed << std::endl;
    file[f] << "saves:                     " << base_settings.number_of_saves << std::endl;
    file[f] << "Speciations per patch:     " << base_settings.speciations_per_patch << std::endl;
    file[f] << "time between saves:        " << base_settings.save_interval() << std::endl;
    // if(check)
    //   file[f] << "bildup_time: " << bildup_time << std::endl;
    // if(check && (bildup_time >= runtime))
    //   file[f] << "WARNING: bildup_time >= runtime" << std::endl;
    file[f] << std::endl;
    file[f] << "============" << std::endl;
    file[f] << std::endl;
    file[f] << "Number of habitats (n²):   " << settings.number_of_habitats() << std::endl;
    file[f] << "Grid length (n):           " << settings.grid_length << std::endl;
    file[f] << "Dispersal range (l):       " << settings.dispersal_range << std::endl;
    if (settings.periodic_boundary_conditions)
        file[f] << "Boundary conditions:       periodic" << std::endl;
    else
        file[f] << "Boundary conditions:       open " << std::endl;
    file[f] << std::endl;
    file[f] << "============" << std::endl;
    file[f] << std::endl;
    file[f] << "Initial dispersal rate:    " << settings.initial_dispersal_rate << std::endl;
    file[f] << "Dispersel variance:        " << settings.dispersal_variance << std::endl;
    file[f] << "Trade-off parameter:       " << settings.zero_crossing << std::endl;
    file[f] << std::endl;
    file[f] << "============" << std::endl;
    file[f] << std::endl;
    file[f] << "Minimum feeding range:     " << settings.min_feeding_range << std::endl;
    file[f] << "Maximum feeding range:     " << settings.max_feeding_range << std::endl;
    file[f] << std::endl;
    file[f] << "Base gain (E0): " << settings.base_gain << std::endl;
    file[f] << "Competition parameter (ξ): " << settings.competition_parameter << std::endl;
    file[f] << "Predation parameter (p):   " << settings.predation_parameter << std::endl;
    file[f] << "Trophic loss (x):          " << settings.trophic_loss << std::endl;
    file[f] << std::endl;
    file[f] << "============" << std::endl;
    file[f] << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_line_habitat_species(resfile_type f,
                                        double time,
                                        size_t habitat,
                                        uint64_t universal_id,
                                        double fitness)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    file[f] << time << "\t" << habitat << "\t" << universal_id << "\t" << fitness << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_line_living_species(resfile_type f,
                                       double time,
                                       uint64_t universal_id,
                                       double first_occurence,
                                       double bodymass,
                                       double feeding_center,
                                       double feeding_range,
                                       double dispersal_rate,
                                       double predator_strength,
                                       size_t population_count,
                                       double mean_trophic_level)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    file[f] << time << "\t" << universal_id << "\t" << first_occurence << "\t" << bodymass << "\t" << feeding_center << "\t" << feeding_range << "\t" << dispersal_rate << "\t" << predator_strength << "\t" << population_count << "\t" << mean_trophic_level << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_line_trophic_levels(resfile_type f,
                                       double time,
                                       std::string habitat,
                                       double dimension,
                                       double mean_trophic_level,
                                       double max_trophic_level)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    file[f] << time << "\t" << habitat << "\t" << dimension << "\t" << mean_trophic_level << "\t" << max_trophic_level << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_line_global_dispersal_info(resfile_type f,
                                              double time,
                                              size_t number_of_species,
                                              size_t number_of_populations,
                                              size_t min_distribution,
                                              double mean_distribution,
                                              size_t max_distribution,
                                              double min_dispersal_rate,
                                              double mean_dispersal_rate_species,
                                              double mean_dispersal_rate_populations,
                                              double max_dispersal_rate,
                                              double min_predator_strength,
                                              double mean_predator_strength_species,
                                              double mean_predator_strength_populations,
                                              double max_predator_strength)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    file[f] << time << "\t" << number_of_species << "\t" << number_of_populations << "\t" << min_distribution << "\t" << mean_distribution << "\t" << max_distribution << "\t" << min_dispersal_rate << "\t" << mean_dispersal_rate_species << "\t" << mean_dispersal_rate_populations << "\t" << max_dispersal_rate << "\t" << min_predator_strength << "\t" << mean_predator_strength_species << "\t" << mean_predator_strength_populations << "\t" << max_predator_strength << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_line_alive_foodwebs(resfile_type f,
                                       double time,
                                       size_t number_of_alive_foodwebs,
                                       double fraction_of_alive_foodwebs)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    file[f] << time << "\t" << number_of_alive_foodwebs << "\t" << fraction_of_alive_foodwebs << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_lifetime_distribution(resfile_type f, double time)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    for (size_t j = 0; j < tl_classes + 1; j++)
    {
        if (lifetime_bins_size[j] > 0)
        {
            file[f] << time << "\t" << get_tl_interval_from_tl_class(j) << "\t" << j << std::endl;

            // take care of the case that the lifetime_bins_size[j] is 0
            for (size_t i = 0; i + 1 < lifetime_bins_size[j]; i++) // -1 prevents a tab at the end
            {
                file[f] << lifetime_bins[j][i] << "\t";
            }

            file[f] << lifetime_bins[j][lifetime_bins_size[j] - 1]; // -1 prevents a tab at the end

            file[f] << std::endl;
            file[f] << std::endl;
        }
    }

    file[f] << std::endl;

    if (!opend)
        file[f].close();
}

void Output::print_LTD_slope(resfile_type f, double time)
{
    if (muted(f))
        return;

    bool opend = file[f].is_open();

    if (!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

    for (size_t tl = 0; tl < tl_classes + 1; tl++)
    {
        file[f] << time << "\t" << get_tl_interval_from_tl_class(tl) << "\t" << calc_LTD_slope(tl) << std::endl;
    }

    file[f] << std::endl;

    if (!opend)
        file[f].close();
}