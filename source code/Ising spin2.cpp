#include "Ising system.hpp"
#include <iostream>
#include <vector>
int main()
{   long long list[3];
    list[0] = 7;
    list[1] = 77;
    list[2] = 777;
    Ising_system mysystem = Ising_system(10, -1.0);

    for (int i = 0; i < 3; i++)
    {
        mysystem.set_state_by_code(list[i]);
        std::vector<int> a=mysystem._state();
        for(int j=0;j<10;j++){
            std::cout<<a[j]<<" ";
        }
        std::cout << mysystem.eval_energy() << std::endl;
        std::cout << mysystem.eval_mz() << std::endl;
    }
    return 0;
}