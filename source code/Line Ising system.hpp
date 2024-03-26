#include "Ising system.hpp"
#include <vector>
#include <cmath>
#include <iostream>

class Line_system:public Ising_system
{   private:
        std::vector<int>system_size;
    public:
        Line_system(const int system_size_spec,const double J_spec):Ising_system(system_size_spec,J_spec),system_size({system_size_spec})
        {
            set_dim(1);
            for(int i=0;i<_n_spins();i++){
                std::vector<int> pos=idx2cod(i);
                set_position(i,pos);
                std::vector<int> NN;
                NN.resize(2);
                std::vector<int>left=pos,right=pos;
                left[0]=(left[0]-1+system_size[0])%system_size[0];
                right[0]=(right[0]+1)%system_size[0];
                NN[0]=cod2idx(left);
                NN[1]=cod2idx(right);
                set_NN(i,NN);
            }
        };
        ~Line_system(){};
        std::vector<int> idx2cod(int idx) override{
            std::vector<int> a={idx};
            return a;
        };
        int cod2idx(std::vector<int> cod) override{
            return cod[0];
        };



};