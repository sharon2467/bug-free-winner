#include "Ising system.hpp"
#include <iostream>
void main()
{
    Ising_system mysystem = Ising_system(10, -1.0);
    int list[3];
    list[1] = 7;
    list[2] = 77;
    list[3] = 777;
    for (int i = 0; i < sizeof(list); i++)
    {
        mysystem.set_state_by_code(7);
        std::cout << mysystem.eval_energy_1D() << std::endl;
        std::cout << mysystem.eval_mz() << std::endl;
    }
}