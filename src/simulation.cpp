#include "simulation.h"

Simulation::Simulation() : m_result(-1)
{

}

Simulation Simulation::create_new(SimulationSettings settings)
{
    Simulation simulation;
    return simulation;
}

Simulation Simulation::load_from_file()
{
    return Simulation();
}

void Simulation::run()
{
    while (m_result == 0) {

    }
}

Simulation::Simulation(std::string path)
{

}