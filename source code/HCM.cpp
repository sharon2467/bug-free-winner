#include "Square Ising system.hpp"
#include <fstream>
void write_in_file(std::ofstream *outputfile, std::vector<double> a)
{
    for (int i = 0; i < size(a); i++)
    {
        *outputfile << a[i] << " ";
    }
    *outputfile << std::endl;
}
int main()
{   std::ofstream outputFile("../data.txt");
    std::ofstream *f = &outputFile;
    std::vector<double> Tlist = {};
    for (double i = 0.05; i < 4; i += 0.05)
    {
        Tlist.push_back(i);
    }
    write_in_file(f, Tlist);
    for (int i = 1; i < 6; i++)
    {

        Square_system a = Square_system(std::vector<int>{i, i}, -1);
        std::vector<double> e , e2,m2 ;
        std::vector<std::vector<double>> result=a.intergrated();
        e=result[1];
        e2=result[2];
        m2=result[0];
        std::vector<double> c(size(Tlist), 0);
        for (int j = 0; j < size(Tlist); j++)
        {
            c[j] = (e2[j]*a._n_spins() - pow(e[j]*a._n_spins(), 2)) / pow(Tlist[j], 2)/a._n_spins();
        }

        // 检查文件是否成功打开
        if (outputFile.is_open())
        {

            write_in_file(f, m2);
            write_in_file(f, c);
            write_in_file(f, e);
        }
        else
        {
            std::cout << "无法打开文件" << std::endl;
        }
    }
}