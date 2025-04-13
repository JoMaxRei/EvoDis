#include "output.h"
#include "easylogging++.h"
#include <cmath>

Output::Output(std::string path, size_t _number_of_habitats, double initial_dispersal_rate, double zero_crossing)
{
    create_file_names(path);

    number_of_habitats = _number_of_habitats;

    tl_classes = 0;

    lifetime_bins_size.push_back(0);

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
        "steps" // 1
        ,
        "living" // 2
        ,
        "TL" // 3
        ,
        "TL_all" // 4
        ,
        "global_info" // 5
        ,
        "alive" // 6
        ,
        "lifetime" // 7
        ,
        "LTD_slope" // 8
        ,
        "abort" // 9
    };

    for (const auto &suffix : suffixes)
    {
        names.push_back(path + "/mc_" + suffix + ".out");
    }
}

bool Output::create_new_files()
{
    bool no_error = true;
    for (size_t i = 0; i < OUT_FILE_COUNT; i++)
    {
        file[i].open(names[i]);
        if (!file[i].is_open())
        {
            LOG(ERROR) << "Could not open file: " << names[i];
            no_error = false;
        }
        file[i].close();
    }

    return no_error;
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

// void Output::get_curr_lifteime_bin(double lifetime)

void Output::update_bins(double lifetime, size_t tl_class)
{
    // Not necessary to update the bins if the output is muted
    if (muted(OUT_LIFETIME) && muted(OUT_LTD_SLOPE))
        return;

    // New species has bigger trophic level than any before
    while (tl_class > tl_classes)
    {

        tl_classes++;

        lifetime_bins_size.push_back(0);

        std::vector<size_t> vec;
        lifetime_bins.push_back(vec);
    }

    double curr_lifetime_bin_d = 0.0;
    size_t curr_lifetime_bin = 0;

    if (lifetime > 0.0)
    {
        curr_lifetime_bin_d = (floor_step(std::log10(lifetime), 1. / static_cast<size_t>(inverted_binsize)) + smallest_lifetime_exponent) * static_cast<size_t>(inverted_binsize);
        if (curr_lifetime_bin_d < 0.0)
            curr_lifetime_bin = 0;
        else
            curr_lifetime_bin = static_cast<size_t>(curr_lifetime_bin_d + machine_epsilon);
    }
    else
        curr_lifetime_bin = 0;

    if (curr_lifetime_bin < lifetime_bins_size[0])
    {
        lifetime_bins[0][curr_lifetime_bin]++;
    }
    else
    {
        for (size_t i = 0; i < curr_lifetime_bin - lifetime_bins_size[0]; i++)
            lifetime_bins[0].push_back(0);

        lifetime_bins[0].push_back(1);
        lifetime_bins_size[0] = curr_lifetime_bin + 1;
    }

    if (tl_class < 1) // Ressource
        return;

    if (curr_lifetime_bin < lifetime_bins_size[tl_class])
    {
        lifetime_bins[tl_class][curr_lifetime_bin]++;
    }
    else
    {
        for (int i = 0; i < curr_lifetime_bin - lifetime_bins_size[tl_class]; i++)
            lifetime_bins[tl_class].push_back(0);

        lifetime_bins[tl_class].push_back(1);
        lifetime_bins_size[tl_class] = curr_lifetime_bin + 1;
    }
}

void Output::print_settings(resfile_type f, BaseSettings base_settings, SimulationSettings settings, bool new_simulation)
{
    if (muted(f))
        return;
    
    bool opend = file[f].is_open();

    if(!opend)
        file[f].open(names[f].c_str(), std::ios::out | std::ios::app);

        file[f] << "===============================================" << std::endl;
        if(new_simulation)
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
        if(settings.periodic_boundary_conditions)
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

    
    if(!opend)
        file[f].close();
}