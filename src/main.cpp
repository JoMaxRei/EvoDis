#include <stdlib.h>

#include "simulation.h"
#include "simulation_settings.h"

int main(int argc, char** argv) {
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