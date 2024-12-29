#include "simulation.h"


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

Simulation::Simulation()
{

}