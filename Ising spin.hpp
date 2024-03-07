class valueException{};
class Ising_spin{
    private:
        int spin;
    public:
        Ising_spin(const int spin_spec=-1):spin(spin_spec){};
        void setvalue(const int value){
            if(value!=1||value!=-1)
                throw valueException();
            else
                spin=value;
        }
        int readvalue() const {
            return spin;
        };
        void set_up(){spin=1;};
        void set_down(){spin=-1;};
        void flip(){spin*=-1;};
        


};
