#include <stdlib.h>

#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

#include "simulation.h"
#include "simulation_settings.h"

int main(int argc, char** argv) {
    // Load config file from disk
    el::Configurations conf("./log.conf");
    // reconfigure loggers to use config
    el::Loggers::reconfigureAllLoggers(conf);
    // log
    LOG(INFO) << "My first info log using default logger";
    LOG(DEBUG) << "My first debug log using default logger";
    LOG(ERROR) << "My first error log using default logger";
    LOG(WARNING) << "My first Warning log using default logger";
    // collect parameters from argc and argv 
    SimulationSettings settings = SimulationSettings::DEFAULT();
    Simulation sim;
    bool load_from_file = false;
    if (load_from_file)
    {
        sim = Simulation::load_from_file(settings);
    }
    else
    {
        sim = Simulation::create_new(settings);
    }
    sim.run();
    return sim.m_result;
}