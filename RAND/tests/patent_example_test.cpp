#include <gtest/gtest.h>
#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include "thermo_backend.hpp"
#include "element_matrix_builder.hpp"
#include <cmath>
#include <iostream>
#include <chrono>

using namespace randflash;
using namespace thermo;

// Test fixture for patent examples
class PatentExampleTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;

    void SetUp() override {
        linSolver = ls::createEigenSolver();
    }

    // Helper: verify element conservation
    bool verifyElementConservation(
        const std::vector<std::vector<double>>& elementMatrix,
        const MultiFlashResult& result,
        const std::vector<double>& expectedElementMoles,
        double tolerance = 1e-6)
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
};

// Patent Example 2: Hydrocarbon cracking multi-component system
TEST_F(PatentExampleTest, Example2_HydrocarbonCracking) {
    std::cout << "\n=== Patent Example 2: Hydrocarbon Cracking Multi-Component System ===" << std::endl;

    std::vector<SpeciesFormula> species = {
        parseFormula("C1", "CH4"),
        parseFormula("C2", "C2H6"),
        parseFormula("C2_1", "C2H4"),
        parseFormula("C3", "C3H8"),
        parseFormula("C3_1", "C3H6"),
        parseFormula("H2", "H2")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    std::vector<double> feedMoles = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "Feed: C3H8 = 1.0 mol" << std::endl;
    std::cout << "Elements: C = " << elementMoles[0] << " mol, H = " << elementMoles[1] << " mol" << std::endl;

    auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C2_1,C3,C3_1,H2", "PR");
    RandFlash flash(*backend, *linSolver);

    FlashInput input;
    input.temperature = 900.0;
    input.pressure = 2.0e5;
    input.feedMoles = feedMoles;

    auto start = std::chrono::high_resolution_clock::now();
    auto result = flash.solveReactive(input, A, 1, 50, 1e-8);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\n--- Results ---" << std::endl;
    std::cout << "Converged: " << (result.success ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Time: " << duration.count() << " ms" << std::endl;

    if (result.success) {
        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles));

        std::cout << "\nFinal composition:" << std::endl;
        const auto& phase = result.phases[0];
        double totalMoles = 0.0;
        for (double n : phase.state.moleNumbers) totalMoles += n;

        std::vector<std::string> names = {"CH4", "C2H6", "C2H4", "C3H8", "C3H6", "H2"};
        for (size_t i = 0; i < phase.state.moleNumbers.size(); ++i) {
            double n = phase.state.moleNumbers[i];
            std::cout << "  " << names[i] << ": " << n << " mol (" 
                      << (n/totalMoles*100.0) << "%)" << std::endl;
        }

        flash.printResult(result);
    }
}

// Patent Example 3: Two-phase reactive system
TEST_F(PatentExampleTest, Example3_TwoPhaseReactive) {
    std::cout << "\n=== Patent Example 3: Two-Phase Reactive System ===" << std::endl;

    std::vector<SpeciesFormula> species = {
        parseFormula("C2", "C2H6"),
        parseFormula("C2_1", "C2H4"),
        parseFormula("H2", "H2")
    };

    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    std::vector<double> feedMoles = {1.0, 0.0, 0.0};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "Feed: C2H6 = 1.0 mol" << std::endl;
    std::cout << "Elements: C = " << elementMoles[0] << " mol, H = " << elementMoles[1] << " mol" << std::endl;

    auto backend = std::make_unique<ThermoPackBackend>("C2,C2_1,H2", "PR");
    RandFlash flash(*backend, *linSolver);

    FlashInput input;
    input.temperature = 300.0;
    input.pressure = 10.0e5;
    input.feedMoles = feedMoles;

    auto start = std::chrono::high_resolution_clock::now();
    auto result = flash.solveReactive(input, A, 2, 50, 1e-8);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\n--- Results ---" << std::endl;
    std::cout << "Converged: " << (result.success ? "Yes" : "No") << std::endl;
    std::cout << "Iterations: " << result.iterations << std::endl;
    std::cout << "Time: " << duration.count() << " ms" << std::endl;
    std::cout << "Phases: " << result.phases.size() << std::endl;

    if (result.success) {
        EXPECT_TRUE(verifyElementConservation(A, result, elementMoles));

        std::vector<std::string> names = {"C2H6", "C2H4", "H2"};
        for (size_t p = 0; p < result.phases.size(); ++p) {
            const auto& phase = result.phases[p];
            std::cout << "\nPhase " << (p+1) << ":" << std::endl;

            double totalMoles = 0.0;
            for (double n : phase.state.moleNumbers) totalMoles += n;

            for (size_t i = 0; i < phase.state.moleNumbers.size(); ++i) {
                double n = phase.state.moleNumbers[i];
                std::cout << "  " << names[i] << ": " << n << " mol (" 
                          << (n/totalMoles*100.0) << "%)" << std::endl;
            }
        }

        flash.printResult(result);
    }
}

int main(int argc, char** argv) {
    randflash::setLogLevel(spdlog::level::debug);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
