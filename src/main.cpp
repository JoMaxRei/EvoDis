#include <stdlib.h>
#include <iostream>

#include "CLI11.hpp"

#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "simulation.h"
#include "simulation_settings.h"

int main(int argc, char **argv)
{
    CLI::App app{"EvoDis - Simulation of a model with evolving dispersal rates"};
    argv = app.ensure_utf8(argv);
    // create subcommands
    CLI::App *create_new = app.add_subcommand("new", "creates a new simulation");
    CLI::App *load_from_file = app.add_subcommand("load", "loads a simulation from file");
    app.require_subcommand(1);

    BaseSettings base_settings = BaseSettings::DEFAULT();
    // define general parameters
    app.add_option("-o,--output", base_settings.output_path, "output path of the simulation")->required();
    app.add_option("-t,--spp", base_settings.speciations_per_patch, "how many events should be simulated on average on each patch before simulation stops")->required();
    app.add_option("-S,--saves", base_settings.number_of_saves, "how many saves during whole simulation")->required();
    app.add_option("-m,--mute", base_settings.muted_outputs, "Outputs to mute (comma-separated): settings,habitat,living,tl,global,alive,ltd,slope,abort")->delimiter(',');

    // define parameters for new sim
    SimulationSettings settings = SimulationSettings::DEFAULT();
    create_new->add_option("-s,--seed", settings.seed, "seed of the simulation")->required();
    create_new->add_option("-g,--gridl", settings.grid_length, "grid length of the simulation; Default: 5");

    create_new->add_option("-d,--drange", settings.dispersal_range, "dispersal range of a species; Default: 1");

    create_new->add_option("-e,--extr", settings.extinction_rate, "random extinction rate; Default: 0");

    create_new->add_option("-i,--inidr", settings.initial_dispersal_rate, "initial dispersal rate; Default: 10.0");
    create_new->add_option("-z,--zeroc", settings.zero_crossing, "dispersal rate at which the predation strength is zero; Default: 1e5");

    // define parameters for old sim
    std::string input_path;
    load_from_file->add_option("-i,--input", input_path, "input path for simulation")->required();

    // parse commands
    CLI11_PARSE(app, argc, argv);

    // Logging configuration
    // Load config file from disk
    el::Configurations conf("./log.conf");
    // reconfigure loggers to use config
    el::Loggers::reconfigureAllLoggers(conf);
    
    // Actual simulation
    Simulation sim;
    if (load_from_file->parsed())
    {
        sim = Simulation(base_settings, input_path);
    }
    else if (create_new->parsed())
    {
        sim = Simulation(base_settings, settings);
    }
    sim.run();
    return sim.m_result;
}