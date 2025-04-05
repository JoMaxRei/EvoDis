#include <stdlib.h>
#include <iostream>

#include "CLI11.hpp"

#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "simulation.h"
#include "simulation_settings.h"

int main(int argc, char **argv)
{
    CLI::App app{"EvoDis - Simulation to prove the developed model"};
    argv = app.ensure_utf8(argv);
    // create subcommands
    CLI::App *create_new = app.add_subcommand("new", "creates a new simulation");
    CLI::App *load_from_file = app.add_subcommand("load", "loads a simulation from file");
    app.require_subcommand(1);

    BaseSettings base_settings;
    // define general parameters
    app.add_option("-o,--output", base_settings.output_path, "output path of the simulation")->required();
    app.add_option("-S,--spp", base_settings.speciations_per_patch, "how many events should be simulated on average on each patch before simulation stops")->required();

    // define parameters for new sim
    SimulationSettings settings = SimulationSettings::DEFAULT();
    create_new->add_option("-s,--seed", settings.seed, "seed of the simulation")->required();
    create_new->add_option("-g,--gridL", settings.grid_length, "grid length of the simulation; Default: 5");

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
    // log
    // LOG(INFO) << "My first info log using default logger";
    // LOG(DEBUG) << "My first debug log using default logger";
    // LOG(ERROR) << "My first error log using default logger";
    // LOG(WARNING) << "My first Warning log using default logger";

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