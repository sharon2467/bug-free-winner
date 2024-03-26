#include <iostream>
#include <vector>
class valueException
{
};
class Ising_spin
{
private:

    int spin;
    std::vector<int> position;
    std::vector<int> NN;
public:
    Ising_spin(const int spin_spec = -1) : spin(spin_spec){};
    void setvalue(const int value)
    {
        if (value != 1 && value != -1)
            throw valueException();
        else
            spin = value;
    }
    int readvalue() const
    {
        return spin;
    };
    void set_up() { spin = 1; };
    void set_down() { spin = -1; };
    void flip() { spin *= -1; };
    void set_dim(int dim)
    {
        position.assign(dim, 0);
    }
    std::vector<int> _position() { return position; };
    std::vector<int> _NN() { return NN; };
    int _NN(const int bondidx) {if(bondidx<NN.size()) return NN[bondidx]; };
    void set_position(std::vector<int> a) { position = a; };
    void set_NN(std::vector<int> a) { NN = a; };
};