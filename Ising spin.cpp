#include "Ising spin.hpp"
#include <iostream>
int main(){
    Ising_spin spin1;
    spin1.create();
    std::cout<<spin1.readvalue()<<std::endl;

}