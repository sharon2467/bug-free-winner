#include "Ising system.hpp"
#include <vector>
#include <cmath>
#include <iostream>
// system_size第一个是行数，第二个是列数，计数从左到右换行,坐标也是{行数，列数}，从0开始,NN是上下左右
class Square_system : public Ising_system
{
private:
    std::vector<int> system_size;

public:
    Square_system(const std::vector<int> system_size_spec, const double J_spec) : Ising_system(system_size_spec[0] * system_size_spec[1], J_spec),
                                                                                  system_size(system_size_spec)
    {
        set_dim(2);
        for (int i = 0; i < _n_spins(); i++)
        {
            std::vector<int> pos = idx2cod(i);
            set_position(i, pos);
            std::vector<int> NN;
            NN.resize(4);
            std::vector<int> up = pos, down = pos, left = pos, right = pos;
            up[0] = (up[0] + 1) % system_size[0];
            down[0] = (down[0] - 1 + system_size[0]) % system_size[0];
            left[1] = (left[1] - 1 + system_size[1]) % system_size[1];
            right[1] = (right[1] + 1) % system_size[1];
            NN[0] = cod2idx(up);
            NN[1] = cod2idx(down);
            NN[2] = cod2idx(left);
            NN[3] = cod2idx(right);
            set_NN(i, NN);
        }
    };
    ~Square_system(){};
    std::vector<int> idx2cod(int idx) override
    {
        std::vector<int> a = {idx / system_size[1], idx % system_size[1]};
        return a;
    };
    int cod2idx(std::vector<int> cod) override
    {
        return cod[0] * system_size[1] + cod[1];
    };
    void print_matrix()
    {
        for (int i = 0; i < n_spins; i++)
        {
            if (i % system_size[1] == 0)
                std::cout << std::endl;
            std::cout << _state(i);
        }
    }
    std::vector<std::vector<int>> _state_matrix()
    {
        std::vector<std::vector<int>> a;
        a.resize(system_size[0]);
        for (int j = 0; j < system_size[0]; j++)
        {
            std::vector<int> b;
            b.resize(system_size[1]);
            for (int i = 0; i < system_size[0]; i++)
            {
                b[i] = spin[j * system_size[0] + i].readvalue();
            }
            a[j] = b;
        }
        return a;
    }
    std::vector<double> mag_sqaure_avg(double Tmin = 0.05, double Tmax = 4, double delta_T = 0.05)
    {
        std::vector<double> Tlist = {};
        for (double i = Tmin; i < Tmax; i += delta_T)
        {
            Tlist.push_back(i);
        }
        std::vector<double> mag_sq_list(size(Tlist), 0), Zlist(size(Tlist), 0);
        int s = _state_code();
        for (int i = 0; i < (maxrep_state + 1) / 2; i++)
        {
            set_state_by_code(i);
            double ground_energy = n_spins * 2 * J;
            for (int j = 0; j < size(mag_sq_list) + 1; j++)
            {
                mag_sq_list[j] += pow(eval_mz(), 2) * exp(-(eval_energy() - ground_energy) / Tlist[j]);
                Zlist[j] += exp(-(eval_energy() - ground_energy) / Tlist[j]);
            }
        }
        for (int j = 0; j < size(mag_sq_list) + 1; j++)
        {
            mag_sq_list[j] /= Zlist[j];
            mag_sq_list[j] /= pow(n_spins, 2);
        }
        set_state_by_code(s);
        return mag_sq_list;
    }
    std::vector<double> energy_pow_avg(int p, double Tmin = 0.05, double Tmax = 4, double delta_T = 0.05)
    {
        std::vector<double> Tlist = {};
        for (double i = Tmin; i < Tmax; i += delta_T)
        {
            Tlist.push_back(i);
        }
        std::vector<double> energy_list(size(Tlist), 0), Zlist(size(Tlist), 0);
        int s = _state_code();
        for (int i = 0; i < (maxrep_state + 1) / 2; i++)
        {
            set_state_by_code(i);
            double ground_energy = n_spins * 2 * J;
            for (int j = 0; j < size(energy_list) + 1; j++)
            {
                energy_list[j] += pow(eval_energy(), p) * exp(-(eval_energy() - ground_energy) / Tlist[j]);
                Zlist[j] += exp(-(eval_energy() - ground_energy) / Tlist[j]);
            }
        }
        for (int j = 0; j < size(energy_list) + 1; j++)
        {
            energy_list[j] /= Zlist[j];
            energy_list[j] /= n_spins;
        }
        set_state_by_code(s);
        return energy_list;
    }
    std::vector<std::vector<double>> intergrated(double Tmin = 0.05, double Tmax = 4, double delta_T = 0.05)
    {
        std::vector<double> Tlist = {};
        for (double i = Tmin; i < Tmax; i += delta_T)
        {
            Tlist.push_back(i);
        }
        std::vector<double> mag_sq_list(size(Tlist), 0), energy_list(size(Tlist), 0), energy_list2(size(Tlist), 0), Zlist(size(Tlist), 0);
        int s = _state_code();
        int energy, mz, previous_energy, previous_mz;
        int ground_energy = n_spins * 2 * J;
        int state;
        for (int i = 0; i < (maxrep_state + 1) / 2; i++)
        {

            if (i == 0)
            {
                set_state_by_code(i);
                energy = eval_energy();
                mz = eval_mz();
            }
            else if (i != 0)
            {
                state = which_state_to_change(i);
                mz = previous_mz - 2 * J * spin[state].readvalue();
                int energy_changed = 0;
                for (int j = 0; j < 4; j++)
                    energy_changed = energy_changed - spin[spin[state]._NN(j)].readvalue() * spin[state].readvalue() * J;
                energy = previous_energy + energy_changed * 2;
                flip_spin(state);
            }                
            previous_energy = energy;
            previous_mz = mz;

            for (int j = 0; j < size(energy_list) + 1; j++)
            {

                energy_list[j] += energy * exp(-(energy - ground_energy) / Tlist[j]);
                mag_sq_list[j] += pow(mz, 2) * exp(-(energy - ground_energy) / Tlist[j]);
                energy_list2[j] += pow(energy, 2) * exp(-(energy - ground_energy) / Tlist[j]);
                Zlist[j] += exp(-(energy - ground_energy) / Tlist[j]);
            }
        }
        for (int j = 0; j < size(energy_list) + 1; j++)
        {
            mag_sq_list[j] /= Zlist[j];
            mag_sq_list[j] /= pow(n_spins, 2);
            energy_list[j] /= Zlist[j];
            energy_list[j] /= n_spins;
            energy_list2[j] /= Zlist[j];
            energy_list2[j] /= n_spins;
        }
        set_state_by_code(s);
        std::vector<std::vector<double>> a;
        a.push_back(mag_sq_list);
        a.push_back(energy_list);
        a.push_back(energy_list2);
        return a;
    }
    int which_state_to_change(int state_idx)
    {
        int i = 0;
        int divisor = 2;
        while (true)
        {
            if (state_idx % divisor == (divisor / 2))
            {
                return i;
            }
            else
            {
                i++;
                divisor *= 2;
            }
        }
    }
};
