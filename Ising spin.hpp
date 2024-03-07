class valueException{};
class Ising_spin{
    private:
        int spin;
    public:
        void create(){
            spin=1;
        }
        void setvalue(const int value){
            if(value!=1||value!=-1)
                throw valueException();
            else
                spin=value;
        }
        int readvalue(){
            return spin;
        }

};
