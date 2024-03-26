#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "Square Ising system.hpp"

TEST_CASE( "Energy Power Average is computed", "[energy_pow_avg]" ) {
    Square_system system=Square_system(std::vector<int>{3,3},-1);

    SECTION("Testing with default parameters") {
        std::vector<double> result = system.energy_pow_avg(2);
        // Add assertions to check the result, e.g., size, values, etc.
        // This will depend on your specific expectations and the behavior of the function
        REQUIRE(result.size() == 80); // Assuming the size should be 80 based on default parameters
    }

    SECTION("Testing with custom parameters") {
        std::vector<double> result = system.energy_pow_avg(3, 0.1, 5, 0.1);
        // Add assertions to check the result
        REQUIRE(result.size() == 50); // Assuming the size should be 50 based on these parameters
    }

    // Add more tests as needed, checking for edge cases, invalid inputs, etc.
}