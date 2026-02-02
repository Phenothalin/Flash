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

    // Use new solveReactive() interface for reactive systems
    auto result = flash.solveReactive(input, A, 1, 50, 1e-8);

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

    // Feed composition (same as TwoPhaseDifferentiatedInit for consistency)
    std::vector<double> feedMoles = {0.5, 0.3, 0.2};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "Initial element moles: C=" << elementMoles[0]
              << ", H=" << elementMoles[1] << std::endl;

    // Setup thermodynamic backend
    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
    RandFlash flash(*backend, *linSolver);

    // Conditions for two-phase region (similar to TwoPhaseDifferentiatedInit)
    FlashInput input;
    input.temperature = 250.0;  // K (two-phase region)
    input.pressure = 2e6;       // Pa (two-phase region)
    input.feedMoles = feedMoles;

    // Use new solveReactive() interface for reactive systems
    auto result = flash.solveReactive(input, A, 2, 50, 1e-8);

    std::cout << "Iterations: " << result.iterations << std::endl;
    EXPECT_TRUE(result.success) << "Flash calculation failed";

    if (result.success) {
        printElementMoles(A, elementNames, result);
        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles))
            << "Element conservation violated";
        flash.printResult(result);
    }
}

// Test 3: Multi-component system with more species
// Verifies element conservation with 5 hydrocarbon species
// Note: PR EOS doesn't model chemical reactions, so composition changes
// are driven by Gibbs energy minimization under element constraints
TEST_F(ReactiveFlashTest, FiveComponentElementConservation) {
    // System: C1, C2, C3, nC4, nC5 (5 hydrocarbons)
    std::vector<SpeciesFormula> species = {
        parseFormula("C1", "CH4"),
        parseFormula("C2", "C2H6"),
        parseFormula("C3", "C3H8"),
        parseFormula("nC4", "C4H10"),
        parseFormula("nC5", "C5H12")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Element matrix:
    //       C1   C2   C3   nC4  nC5
    // C  [  1    2    3    4    5  ]
    // H  [  4    6    8   10   12  ]

    // Mixed initial composition (with light components for two-phase region)
    std::vector<double> feedMoles = {0.3, 0.25, 0.2, 0.15, 0.1};

    // Compute element moles from feed
    auto elementMoles = computeElementMoles(A, feedMoles);
    std::cout << "\n=== Test 3: Five Component Element Conservation ===" << std::endl;
    std::cout << "Initial species moles: [";
    for (size_t i = 0; i < feedMoles.size(); ++i) {
        std::cout << feedMoles[i] << (i < feedMoles.size()-1 ? ", " : "");
    }
    std::cout << "]" << std::endl;
    std::cout << "Initial element moles: C=" << elementMoles[0]
              << ", H=" << elementMoles[1] << std::endl;

    // Setup thermodynamic backend
    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3,nC4,nC5", "PR");
    RandFlash flash(*backend, *linSolver);

    // Prepare flash input (two-phase region conditions)
    FlashInput input;
    input.temperature = 250.0;  // K (two-phase region)
    input.pressure = 2e6;       // Pa (two-phase region)
    input.feedMoles = feedMoles;

    // Use new solveReactive() interface for reactive systems
    auto result = flash.solveReactive(input, A, 2, 50, 1e-8);

    std::cout << "Iterations: " << result.iterations << std::endl;
    EXPECT_TRUE(result.success) << "Flash calculation failed";

    if (result.success) {
        // Print final species moles
        std::cout << "Final species moles: [";
        for (size_t i = 0; i < result.phases[0].state.moleNumbers.size(); ++i) {
            std::cout << result.phases[0].state.moleNumbers[i];
            if (i < result.phases[0].state.moleNumbers.size()-1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;

        printElementMoles(A, elementNames, result);
        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles))
            << "Element conservation violated";

        // Check that composition changed (algorithm did work)
        bool compositionChanged = false;
        for (size_t i = 0; i < feedMoles.size(); ++i) {
            double diff = std::abs(result.phases[0].state.moleNumbers[i] - feedMoles[i]);
            if (diff > 1e-6) {
                compositionChanged = true;
                break;
            }
        }
        EXPECT_TRUE(compositionChanged)
            << "Expected composition to change during equilibration";

        flash.printResult(result);
    }
}

// Test 4: Zero mole component generation (Problem 3 fix verification)
// Verifies that species with zero initial moles can be generated
TEST_F(ReactiveFlashTest, ZeroMoleComponentGeneration) {
    // System: C1, C2, C3, nC4, nC5
    std::vector<SpeciesFormula> species = {
        parseFormula("C1", "CH4"),
        parseFormula("C2", "C2H6"),
        parseFormula("C3", "C3H8"),
        parseFormula("nC4", "C4H10"),
        parseFormula("nC5", "C5H12")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Feed: C1 starts at ZERO (should be generated during equilibration)
    std::vector<double> feedMoles = {0.0, 0.3, 0.2, 0.1, 0.1};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "\n=== Test 4: Zero Mole Component Generation ===" << std::endl;
    std::cout << "Initial C1 moles: " << feedMoles[0] << " (should be zero)" << std::endl;

    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3,nC4,nC5", "PR");
    RandFlash flash(*backend, *linSolver);

    FlashInput input;
    input.temperature = 300.0;
    input.pressure = 1e5;
    input.feedMoles = feedMoles;

    // Use new solveReactive() interface for reactive systems
    auto result = flash.solveReactive(input, A, 1, 50, 1e-8);

    EXPECT_TRUE(result.success) << "Flash calculation failed";

    if (result.success) {
        double final_C1 = result.phases[0].state.moleNumbers[0];
        std::cout << "Final C1 moles: " << final_C1 << std::endl;

        // C1 should be generated (non-zero)
        EXPECT_GT(final_C1, 1e-6) << "C1 should be generated from zero initial value";

        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles))
            << "Element conservation violated";
    }
}

// Test 5: Multi-phase differentiated initialization (Problem 1 fix verification)
// Verifies that two-phase initialization produces different compositions
TEST_F(ReactiveFlashTest, TwoPhaseDifferentiatedInit) {
    std::vector<SpeciesFormula> species = {
        parseFormula("C1", "CH4"),
        parseFormula("C2", "C2H6"),
        parseFormula("C3", "C3H8")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    std::vector<double> feedMoles = {0.5, 0.3, 0.2};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "\n=== Test 5: Two-Phase Differentiated Initialization ===" << std::endl;

    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
    RandFlash flash(*backend, *linSolver);

    FlashInput input;
    input.temperature = 250.0;  // Two-phase region
    input.pressure = 2e6;
    input.feedMoles = feedMoles;

    // Use new solveReactive() interface for reactive systems
    auto result = flash.solveReactive(input, A, 2, 50, 1e-8);

    EXPECT_TRUE(result.success) << "Flash calculation failed";
    EXPECT_EQ(result.numPhases(), 2) << "Should have 2 phases";
    EXPECT_GT(result.iterations, 1) << "Should require more than 1 iteration";

    if (result.success && result.numPhases() == 2) {
        // Check that phases have different compositions
        double composition_diff = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            double x1 = result.phases[0].x[i];
            double x2 = result.phases[1].x[i];
            composition_diff += std::abs(x1 - x2);
        }

        std::cout << "Composition L1 distance: " << composition_diff << std::endl;
        EXPECT_GT(composition_diff, 0.01) << "Phases should have different compositions";

        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles))
            << "Element conservation violated";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
