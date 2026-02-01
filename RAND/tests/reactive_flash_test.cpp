#include <gtest/gtest.h>
#include "rand_flash.hpp"
#include "thermo_backend.hpp"
#include "element_matrix_builder.hpp"
#include <cmath>
#include <iostream>
#include <numeric>

using namespace randflash;
using namespace thermo;

// Test fixture for reactive flash calculations
class ReactiveFlashTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;

    void SetUp() override {
        linSolver = ls::createEigenSolver();
    }

    // Helper: verify element conservation across all phases
    bool verifyElementConservation(
        const std::vector<std::vector<double>>& elementMatrix,
        const MultiFlashResult& result,
        const std::vector<double>& expectedElementMoles,
        double tolerance = 1e-6)
    {
        const size_t E = elementMatrix.size();
        const size_t C = elementMatrix[0].size();

        std::vector<double> computedElementMoles(E, 0.0);

        // Sum over all phases
        for (const auto& phase : result.phases) {
            for (size_t e = 0; e < E; ++e) {
                for (size_t i = 0; i < C; ++i) {
                    computedElementMoles[e] += elementMatrix[e][i] * phase.state.moleNumbers[i];
                }
            }
        }

        // Check each element
        bool success = true;
        for (size_t e = 0; e < E; ++e) {
            double diff = std::abs(computedElementMoles[e] - expectedElementMoles[e]);
            double scale = std::max(std::abs(expectedElementMoles[e]), 1.0);
            if (diff > tolerance * scale) {
                std::cout << "Element " << e << " conservation failed: "
                          << "computed=" << computedElementMoles[e]
                          << ", expected=" << expectedElementMoles[e]
                          << ", diff=" << diff << std::endl;
                success = false;
            }
        }
        return success;
    }

    // Helper: print element moles
    void printElementMoles(
        const std::vector<std::vector<double>>& elementMatrix,
        const std::vector<std::string>& elementNames,
        const MultiFlashResult& result)
    {
        const size_t E = elementMatrix.size();
        const size_t C = elementMatrix[0].size();

        std::vector<double> computedElementMoles(E, 0.0);
        for (const auto& phase : result.phases) {
            for (size_t e = 0; e < E; ++e) {
                for (size_t i = 0; i < C; ++i) {
                    computedElementMoles[e] += elementMatrix[e][i] * phase.state.moleNumbers[i];
                }
            }
        }

        std::cout << "Element moles after equilibrium:" << std::endl;
        for (size_t e = 0; e < E; ++e) {
            std::cout << "  " << elementNames[e] << ": " << computedElementMoles[e] << std::endl;
        }
    }
};

// Test 1: Single-phase reactive system with element matrix
// Key: disable stability test, force single phase
// Verifies element conservation with non-identity element matrix
TEST_F(ReactiveFlashTest, SinglePhaseElementConservation) {
    // System: C1 (CH4), C2 (C2H6), C3 (C3H8)
    std::vector<SpeciesFormula> species = {
        parseFormula("C1", "CH4"),
        parseFormula("C2", "C2H6"),
        parseFormula("C3", "C3H8")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Element matrix:
    //      C1   C2   C3
    // C  [  1    2    3  ]
    // H  [  4    6    8  ]

    // Feed composition (initial species moles)
    std::vector<double> feedMoles = {0.5, 0.3, 0.2};

    // Compute element moles from feed
    auto elementMoles = computeElementMoles(A, feedMoles);
    std::cout << "Initial element moles: C=" << elementMoles[0]
              << ", H=" << elementMoles[1] << std::endl;

    // Setup thermodynamic backend
    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
    RandFlash flash(*backend, *linSolver);

    // Prepare flash input
    FlashInput input;
    input.temperature = 300.0;  // K
    input.pressure = 1e5;       // Pa
    input.feedMoles = feedMoles;

    // IMPORTANT: For reactive systems, disable stability test
    // and manually specify phase count
    SolveOptions options;
    options.enable_stability_test = false;
    options.forced_phase_count = 1;

    auto result = flash.solve(input, A, options, 50, 1e-8);

    std::cout << "Iterations: " << result.iterations << std::endl;
    EXPECT_TRUE(result.success) << "Flash calculation failed";

    if (result.success) {
        printElementMoles(A, elementNames, result);
        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles))
            << "Element conservation violated";
        flash.printResult(result);
    }
}

// Test 2: Two-phase reactive system (VLE with element conservation)
// Verifies element conservation across vapor and liquid phases
TEST_F(ReactiveFlashTest, TwoPhaseElementConservation) {
    // System: C1, C2, C3
    std::vector<SpeciesFormula> species = {
        parseFormula("C1", "CH4"),
        parseFormula("C2", "C2H6"),
        parseFormula("C3", "C3H8")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Feed composition
    std::vector<double> feedMoles = {0.6, 0.3, 0.1};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "Initial element moles: C=" << elementMoles[0]
              << ", H=" << elementMoles[1] << std::endl;

    // Setup thermodynamic backend
    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
    RandFlash flash(*backend, *linSolver);

    // Conditions for two-phase region
    FlashInput input;
    input.temperature = 200.0;  // K (low temp)
    input.pressure = 3e6;       // Pa (high pressure)
    input.feedMoles = feedMoles;

    // Disable stability test, force two phases
    SolveOptions options;
    options.enable_stability_test = false;
    options.forced_phase_count = 2;

    auto result = flash.solve(input, A, options, 50, 1e-8);

    std::cout << "Iterations: " << result.iterations << std::endl;
    EXPECT_TRUE(result.success) << "Flash calculation failed";

    if (result.success) {
        printElementMoles(A, elementNames, result);
        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles))
            << "Element conservation violated";
        flash.printResult(result);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
