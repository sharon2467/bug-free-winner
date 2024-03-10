#include "Ising spin.hpp"
#include <vector>
#include <cmath>
#include <iostream>
class Ising_system
{
private:
    const double J;
    const int n_spins;
    const long long maxrep_state;
    std::vector<Ising_spin> spin;

public:
    Ising_system(const int n_spins_spec, const double J1) : J(J1), n_spins(n_spins_spec),maxrep_state(static_cast<long long>(std::pow(2, n_spins)) - 1)
    {
        spin.resize(n_spins);
    };
    virtual ~Ising_system(){};

    double _J() const { return J; };
    int _n_spins() const { return n_spins; };
    long long _maxrep_state() const { return maxrep_state; };

    int _sz(const int site_idx) const { return spin[site_idx].readvalue(); };
    void set_up_spin(const int site_idx) { spin[site_idx].set_up(); };
    void set_dw_spin(const int site_idx) { spin[site_idx].set_down(); };
    void set_spin(const int site_idx, int s_spec) { spin[site_idx].setvalue(s_spec); };
    void flip_spin(const int site_idx) { spin[site_idx].flip(); };
    void set_state_by_code(long long rep_state)
    {
        int code_rest = rep_state;
        if (rep_state < maxrep_state)
        {
            for (int i = 0; i < n_spins; i++)
            {
                spin[i].set_down();
            }
            int i = 0;
            
            while (code_rest != 0)
            {

                spin[i].setvalue((code_rest % 2 ? 1 : -1));
                code_rest = (code_rest - (code_rest % 2)) / 2;
                i++;
            }
        }
    }
    std::vector<int> return_state() const
    {
        std::vector<int> a;
        a.resize(n_spins);
        for (int i = 0; i < n_spins; i++)
        {
            a[i] = spin[i].readvalue();
        }
        return a;
    }
    long long return_state_code() const
    {
        long long code = 0;
        long long multiplier = 1;
        for (int i = 0; i < n_spins; i++)
        {
            code += (spin[i].readvalue() + 1) * multiplier;
            multiplier *= 2;
        }
        return code / 2;
    }
    double eval_mz() const
    {
        double mz = 0;
        for (int i = 0; i < n_spins; i++)
        {

            mz = mz + J * spin[i].readvalue();
        }
        return mz;
    }
    double eval_energy_1D() const
    {
        double energy = 0;
        energy += J * spin[0].readvalue() * spin[n_spins - 1].readvalue();
        for (int i = 1; i < n_spins; i++)
        {
            energy += J * spin[i].readvalue() * spin[i - 1].readvalue();
        }
        return energy;
    }
};