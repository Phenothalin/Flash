#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include "rand_flash.hpp"
#include "thermo_backend.hpp"
#include "rand_logger.hpp"

class PhaseDeterminationTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;
    std::ofstream dataFile;

    void SetUp() override {
        linSolver = ls::createEigenSolver();
        std::ifstream checkFile("/Users/madao/Desktop/2606/thesis_data/phase_determination.csv");
        bool fileExists = checkFile.good();
        checkFile.close();

        dataFile.open("/Users/madao/Desktop/2606/thesis_data/phase_determination.csv", std::ios::app);
        if (!fileExists) {
            dataFile << "TestCase,System,T_K,P_bar,Phases,Beta,TPD_min,Iterations,Status\n";
        }
    }

    void TearDown() override {
        dataFile.close();
    }

    void RecordResult(const std::string& testName, const std::string& system,
                     double T, double P, const randflash::MultiFlashResult& result, double tpd_min = 0.0) {
        std::string beta_str;
        if (result.numPhases() == 1) {
            beta_str = result.phases[0].state.phaseFlag == 1 ? "1.0" : "0.0";
        } else {
            for (size_t p = 0; p < result.numPhases(); ++p) {
                if (p > 0) beta_str += ";";
                beta_str += std::to_string(result.beta(p));
            }
        }

        dataFile << testName << "," << system << "," << T << "," << P/1e5 << ","
                 << result.numPhases() << "," << beta_str << "," << tpd_min << ","
                 << result.iterations << "," << (result.success ? "CONVERGED" : "FAILED") << "\n";
    }
};

// Case 17: CH₄-C₂H₆ Temperature Sweep
TEST_F(PhaseDeterminationTest, Case17_CH4_C2H6_TempSweep) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2", "SRK");
    randflash::RandFlash flash(*backend, *linSolver);
    std::vector<double> feed = {0.5, 0.5};

    std::vector<std::vector<double>> A(2, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 2; ++i) A[i][i] = 1.0;

    for (double T : {300.0, 280.0, 260.0, 240.0, 230.0}) {
        randflash::FlashInput input{T, 30e5, feed};
        randflash::SolveOptions opts;
        opts.phase_determination_strategy = "stable";

        auto result = flash.solve(input, A, opts);
        flash.printResult(result);
        RecordResult("Case17_TempSweep", "C1-C2", T, 30e5, result);
    }
}

// Case 18: N₂-CH₄ Pressure Sweep
TEST_F(PhaseDeterminationTest, Case18_N2_CH4_PressureSweep) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("N2,C1", "SRK");
    randflash::RandFlash flash(*backend, *linSolver);
    std::vector<double> feed = {0.3, 0.7};

    std::vector<std::vector<double>> A(2, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 2; ++i) A[i][i] = 1.0;

    for (double P_bar : {1.0, 10.0, 30.0, 50.0, 100.0}) {
        randflash::FlashInput input{120.0, P_bar * 1e5, feed};
        randflash::SolveOptions opts;
        opts.phase_determination_strategy = "stable";

        auto result = flash.solve(input, A, opts);
        RecordResult("Case18_PressureSweep", "N2-C1", 120.0, P_bar * 1e5, result);
    }
}

// Case 19: H₂O-nC₆ LLE
TEST_F(PhaseDeterminationTest, Case19_H2O_nC6_LLE) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("H2O,nC6", "SRK");
    randflash::RandFlash flash(*backend, *linSolver);
    std::vector<double> feed = {0.5, 0.5};

    std::vector<std::vector<double>> A(2, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 2; ++i) A[i][i] = 1.0;

    randflash::FlashInput input{298.0, 1e5, feed};
    randflash::SolveOptions opts;
    opts.phase_determination_strategy = "stable";

    auto result = flash.solve(input, A, opts);
    RecordResult("Case19_LLE", "H2O-nC6", 298.0, 1e5, result);
}

// Case 20: H₂O-CH₄-nC₆ VLLE
TEST_F(PhaseDeterminationTest, Case20_H2O_CH4_nC6_VLLE) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("H2O,C1,nC6", "SRK");
    randflash::RandFlash flash(*backend, *linSolver);
    std::vector<double> feed = {0.3, 0.5, 0.2};

    std::vector<std::vector<double>> A(3, std::vector<double>(3, 0.0));
    for (size_t i = 0; i < 3; ++i) A[i][i] = 1.0;

    randflash::FlashInput input{298.0, 50e5, feed};
    randflash::SolveOptions opts;
    opts.phase_determination_strategy = "stable";

    auto result = flash.solve(input, A, opts);
    RecordResult("Case20_VLLE", "H2O-C1-nC6", 298.0, 50e5, result);
}

// Case 21: CO₂-C₃H₈ Near-Critical
TEST_F(PhaseDeterminationTest, Case21_CO2_C3H8_NearCritical) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO2,C3", "SRK");
    randflash::RandFlash flash(*backend, *linSolver);
    std::vector<double> feed = {0.5, 0.5};

    std::vector<std::vector<double>> A(2, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 2; ++i) A[i][i] = 1.0;

    double T_crit = 310.0;
    for (double dT : {-20.0, -10.0, -5.0, -2.0, 0.0}) {
        randflash::FlashInput input{T_crit + dT, 70e5, feed};
        randflash::SolveOptions opts;
        opts.phase_determination_strategy = "stable";

        auto result = flash.solve(input, A, opts);
        RecordResult("Case21_NearCritical", "CO2-C3", T_crit + dT, 70e5, result);
    }
}

// Case 22: CH₄-C₂H₆ Phase Disappearance
TEST_F(PhaseDeterminationTest, Case22_CH4_C2H6_PhaseDisappearance) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2", "SRK");
    randflash::RandFlash flash(*backend, *linSolver);
    std::vector<double> feed = {0.5, 0.5};

    std::vector<std::vector<double>> A(2, std::vector<double>(2, 0.0));
    for (size_t i = 0; i < 2; ++i) A[i][i] = 1.0;

    for (double T : {220.0, 230.0, 240.0, 250.0}) {
        randflash::FlashInput input{T, 30e5, feed};
        randflash::SolveOptions opts;
        opts.phase_determination_strategy = "stable";

        auto result = flash.solve(input, A, opts);
        RecordResult("Case22_PhaseDisappear", "C1-C2", T, 30e5, result);
    }
}

int main(int argc, char **argv) {
    randflash::setLogLevel(spdlog::level::info);
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "PhaseDeterminationTest.Case17_CH4_C2H6_TempSweep";
    return RUN_ALL_TESTS();
}
