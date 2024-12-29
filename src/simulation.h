#ifndef SIMULATION_H_
#define SIMULATION_H_

class Simulation
{
public:
    static Simulation create_new();
    static Simulation load_from_file();
    void run();

    int m_result;
};

#endif