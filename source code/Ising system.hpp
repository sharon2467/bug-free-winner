#include "Ising spin.hpp"
#include <vector>
#include <cmath>
#include <iostream>
class OutOfSystemBoundException
{
};
class Ising_system
{

protected:
    const double J;
    const int n_spins;
    const long long maxrep_state;
    std::vector<Ising_spin> spin;

public:
    Ising_system(const int n_spins_spec, const double J1) : J(J1), n_spins(n_spins_spec), maxrep_state(static_cast<long long>(std::pow(2, n_spins)) - 1)
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
    virtual std::vector<int> idx2cod(int idx) { return std::vector<int>{0}; };
    virtual int cod2idx(std::vector<int>) { return 0; };
    void set_state_by_code(long long rep_state)

    {
        int code_rest = rep_state;
        if (rep_state <= maxrep_state)
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
    std::vector<int> _state() const
    {
        std::vector<int> a;
        a.resize(n_spins);
        for (int i = 0; i < n_spins; i++)
        {
            a[i] = spin[i].readvalue();
        }
        return a;
    }
    int _state(int site_idx) const
    {
        if (site_idx < n_spins)
            return spin[site_idx].readvalue();
        else
            throw OutOfSystemBoundException();
    }
    long long _state_code() const
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
    void set_dim(const int dim)
    {
        for (auto &each : spin)
            each.set_dim(dim);
    };
    std::vector<int> _spin_position(int site_idx)
    {
        if (site_idx < n_spins)
            return spin[site_idx]._position();
        else
            throw OutOfSystemBoundException();
    };
    std::vector<int> _spin_NN(int site_idx)
    {
        if (site_idx < n_spins)
            return spin[site_idx]._NN();
        else
            throw OutOfSystemBoundException();
    };
    int _spin_NN(int site_idx, int bond_idx)
    {
        if (site_idx < n_spins)
            return spin[site_idx]._NN(bond_idx);
        else
            throw OutOfSystemBoundException();
    };
    void set_NN(int site_idx, const std::vector<int> NN)
    {
        if (site_idx < n_spins)
            spin[site_idx].set_NN(NN);
        else
            throw OutOfSystemBoundException();
    }
    void set_position(int site_idx, const std::vector<int> position)
    {
        if (site_idx < n_spins)
            spin[site_idx].set_position(position);
        else
            throw OutOfSystemBoundException();
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
    double eval_energy()
    {
        double energy = 0;
        for (int i = 0; i < n_spins; i++)
        {
            for (int j = 0; j < spin[i]._NN().size(); j++)
                energy = energy + spin[spin[i]._NN(j)].readvalue() * spin[i].readvalue() * J;
        }
        return energy/2;
    }

};