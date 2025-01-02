#include <stdlib.h>

#include "CLI11.hpp"

#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "simulation.h"
#include "simulation_settings.h"

int main(int argc, char** argv) {
    CLI::App app{"EvoDis - Simulation to prove the developed model"};
    argv = app.ensure_utf8(argv);
    // create subcommands
    CLI::App *create_new = app.add_subcommand("new", "creates a new simulation");
    CLI::App *load_from_file = app.add_subcommand("load", "loads a simulation from file");
    app.require_subcommand(1);

    // define general parameters
    std::string output_path;
    app.add_option("-o,--output", output_path, "output path of the simulation")->required();
    uint64_t speciations_per_habitat;
    app.add_option("-sph,--speciations-per-habitat", speciations_per_habitat, "how many events should be simulated before simulation stops")->required();

    // define parameters for new sim
    uint64_t seed;
    create_new->add_option("-s,--seed", seed, "seed of the simulation")->required();

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
    SimulationSettings settings = SimulationSettings::DEFAULT();
    Simulation sim;
    if (load_from_file->parsed()) {
        sim = Simulation::load_from_file(settings);
    }
    else if(create_new->parsed()) {
        sim = Simulation::create_new(settings);
    }
    sim.run();
    return sim.m_result;
}