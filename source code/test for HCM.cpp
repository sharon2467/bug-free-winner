#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "Square Ising system.hpp"

TEST_CASE("Test vector sizes", "[vector_sizes]") {
    std::vector<double> Tlist;
    for (double i = 0.05; i < 4; i += 0.05)
    {
        Tlist.push_back(i);
    }

    REQUIRE(Tlist.size() == 80); // 4/0.05 - 1 = 79

    for (int i = 1; i < 3; i++)
    {
        Square_system a = Square_system(std::vector<int>{i, i}, -1);
        std::vector<double> e = a.energy_pow_avg(1), e2 = a.energy_pow_avg(2),
                            m2 = a.mag_sqaure_avg();
        std::vector<double> c(size(Tlist), 0);

        REQUIRE(e.size() == Tlist.size());
        REQUIRE(e2.size() == Tlist.size());
        REQUIRE(m2.size() == Tlist.size());
        REQUIRE(c.size() == Tlist.size());
    }
}