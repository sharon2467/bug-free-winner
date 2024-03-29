#define CATCH_CONFIG_MAIN  
#include "catch.hpp"
#include "Square Ising system.hpp"
#include <vector>

TEST_CASE("Ising System Tests", "[IsingSystem]") {
    long long list[3] = {7, 77, 777};
    std::vector<int>s={2,5};
    Square_system mysystem(s, -1.0);

    SECTION("Test State and Energy Evaluation") {
        std::vector<std::vector<int>> expectedStates = {{1,1,1,-1,-1,-1,-1,-1,-1,-1}, {1,-1,1,1,-1,-1,1,-1,-1,-1}, {1,-1,-1,1,-1,-1,-1,-1,1,1}};
        std::vector<double> expectedEnergies = {-4,8,0};
        std::vector<int> expectedMz = {4,2,2};

        for (int i = 0; i < 3; i++) {
            mysystem.set_state_by_code(list[i]);
            std::vector<int> actualState = mysystem._state();
            std::vector<std::vector<int>> matrix=mysystem._state_matrix();
            mysystem.print_matrix();
            // 检查状态是否符合预期
            REQUIRE(actualState == expectedStates[i]);

            // 检查能量评估是否符合预期
            double actualEnergy = mysystem.eval_energy();
            REQUIRE(actualEnergy == Approx(expectedEnergies[i]));

            // 检查mz值是否符合预期
            int actualMz = mysystem.eval_mz();
            REQUIRE(actualMz == expectedMz[i]);
        }
    }
}
