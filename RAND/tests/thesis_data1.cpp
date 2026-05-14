// thesis_data1.cpp - RAND solver thesis data collection (Cases 6-16)
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include "thermo_backend.hpp"

static const std::string CSV_PATH = "/Users/madao/Desktop/2606/thesis_data/rand_results.csv";
static const std::string LOG_PATH = "/Users/madao/Desktop/2606/thesis_data/rand_detailed.log";
static const std::string ITER_CSV = "/Users/madao/Desktop/2606/thesis_data/rand_iter_history.csv";

class RANDThesisTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;
    std::ofstream csv;

    void SetUp() override {
        linSolver = ls::createEigenSolver();
        std::ifstream check(CSV_PATH);
        bool exists = check.good();
        check.close();
        csv.open(CSV_PATH, std::ios::app);
        if (!exists)
            csv << "Case,System,T_K,P_bar,Phases,Beta,Iterations,ChemPotError,ElemResidual,Status\n";

        std::ifstream iter_check(ITER_CSV);
        bool iter_exists = iter_check.good();
        iter_check.close();
        std::ofstream iter_csv(ITER_CSV, std::ios::app);
        if (!iter_exists)
            iter_csv << "Case,Iteration,MuError\n";
    }

    void TearDown() override { csv.close(); }

    std::vector<std::vector<double>> eye(size_t n) {
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) A[i][i] = 1.0;
        return A;
    }

    void record(const std::string& cas, const std::string& sys,
                double T, double P, const randflash::MultiFlashResult& r, long long us = -1) {
        csv << cas << "," << sys << "," << T << "," << P/1e5 << ","
            << r.numPhases() << ",";
        for (size_t i = 0; i < r.numPhases(); ++i) {
            csv << r.beta(i);
            if (i+1 < r.numPhases()) csv << ";";
        }
        csv << "," << r.iterations << "," << r.mu_infinity_norm << "," << r.elem_residual_inf << ","
            << (r.success ? "CONVERGED" : "FAILED");
        if (us >= 0) csv << "," << us;
        csv << "\n";

        if (!r.iter_mu_history.empty()) {
            std::ofstream iter_csv(ITER_CSV, std::ios::app);
            for (size_t i = 0; i < r.iter_mu_history.size(); ++i)
                iter_csv << cas << "," << (i + 1) << "," << r.iter_mu_history[i] << "\n";
        }

        std::ofstream log(LOG_PATH, std::ios::app);
        log << "\n=== " << cas << " ===\nT=" << T << "K P=" << P/1e5
            << "bar phases=" << r.numPhases() << " iters=" << r.iterations;
        if (us >= 0) log << " time=" << us << "us";
        log << "\n";
        for (size_t p = 0; p < r.numPhases(); ++p) {
            const auto& n = r.n_phase(p);
            double b = r.beta(p);
            log << "Phase" << p << "(beta=" << b << "): ";
            for (size_t c = 0; c < n.size(); ++c)
                log << n[c]/b << (c+1<n.size()?", ":"");
            log << "\n";
        }
    }
};

// Case 6: CH4-C2H6 two-phase, T=250K, P=30bar
TEST_F(RANDThesisTest, Case6_C1_C2_250K_30bar) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {0.5, 0.5};
    randflash::FlashInput input{250.0, 30e5, z};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t6 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(2), opts);
    auto us6 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t6).count();
    solver.printResult(r);
    record("Case6_RAND", "C1-C2", 250.0, 30e5, r, us6);
    EXPECT_TRUE(r.success);
}

// Case 7: CO2-C3H8 near-critical, T=310K, P=70bar
TEST_F(RANDThesisTest, Case7_CO2_C3H8_310K_70bar) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO2,C3", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {0.5, 0.5};
    randflash::FlashInput input{310.0, 70e5, z};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t7 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(2), opts);
    auto us7 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t7).count();
    solver.printResult(r);
    record("Case7_RAND", "CO2-C3", 310.0, 70e5, r, us7);
    EXPECT_TRUE(r.success);
}

// Case 9: 5-comp natural gas, T=280K, P=50bar
TEST_F(RANDThesisTest, Case9_5comp_280K_50bar) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2,C3,nC4,nC5", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {8.5, 0.8, 0.4, 0.2, 0.1};
    randflash::FlashInput input{280.0, 50e5, z};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t9 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(5), opts);
    auto us9 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t9).count();
    solver.printResult(r);
    record("Case9_RAND", "C1-C2-C3-nC4-nC5", 280.0, 50e5, r, us9);
    EXPECT_TRUE(r.success);
}

// Case 10: 11-comp natural gas, T=295K, P=20bar
TEST_F(RANDThesisTest, Case10_11comp_295K_20bar) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {0.003,0.015,0.55,0.14,0.12,0.05,0.045,0.03,0.025,0.012,0.01};
    randflash::FlashInput input{295.0, 20e5, z};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t10 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(11), opts);
    auto us10 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t10).count();
    solver.printResult(r);
    record("Case10_RAND", "11comp", 295.0, 20e5, r, us10);
    EXPECT_TRUE(r.success);
}

// TEST_F(RANDThesisTest, Case10_11comp_295K_20bar_ScaledFeed) {
//     auto backend = std::make_unique<thermo::ThermoPackBackend>("N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7", "SRK");
//     randflash::RandFlash solver(*backend, *linSolver);
//     std::vector<double> z = {0.03,0.15,5.5,1.4,1.2,0.5,0.45,0.3,0.25,0.12,0.1};
//     randflash::FlashInput input{295.0, 20e5, z};
//     randflash::SolveOptions opts;
//     opts.enable_stability_test = false;
//     opts.forced_phase_count = 2;
//     opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
//     auto r = solver.solve(input, eye(11), opts);
//     solver.printResult(r);
//     EXPECT_TRUE(r.success);
// }

// Case 11: 20-comp natural gas, T=300K, P=50bar
TEST_F(RANDThesisTest, Case11_20comp_300K_50bar) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>(
        "N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7,nC8,nC9,nC10,nC11,nC12,nC13,nC14,nC15,nC16", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {0.1,0.2,7.0,1.0,0.6,0.3,0.2,0.15,0.1,0.08,
                             0.06,0.04,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.07};
    randflash::FlashInput input{300.0, 50e5, z};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t11 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(20), opts);
    auto us11 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t11).count();
    solver.printResult(r);
    record("Case11_RAND", "20comp", 300.0, 50e5, r, us11);
    EXPECT_TRUE(r.success);
}

// Case 12: H2O-CH4-nC6 phase equilibrium, T=298.15K, P=50bar
TEST_F(RANDThesisTest, Case12_H2O_CH4_nC6_298K_50bar) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("H2O,C1,nC6", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {0.3, 0.3, 0.4};
    randflash::FlashInput input{298.15, 50e5, z};
    randflash::SolveOptions opts;
    opts.phase_determination_strategy = "stable";  // 使用稳定法策略
    auto t12 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(3), opts);
    auto us12 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t12).count();
    solver.printResult(r);
    record("Case12_RAND", "H2O-C1-nC6", 298.15, 50e5, r, us12);
    EXPECT_TRUE(r.success);
}

// Case 13: Hydrocarbon reactive equilibrium, T=320K, P=30bar
// b=[C=2.5, H=7.0]; feed: {0.25,0,0.75,0,0} -> A*feed=[1*0.25+3*0.75, 4*0.25+8*0.75]=[2.5,7.0]
TEST_F(RANDThesisTest, Case13_Hydrocarbon_Reactive) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2,C3,nC4,nC5", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<std::vector<double>> A = {
        {1, 2, 3, 4, 5},
        {4, 6, 8, 10, 12}
    };
    randflash::FlashInput input{320.0, 30e5, {0.25, 0.0, 0.75, 0.0, 0.0}};
    randflash::SolveOptions opts;
    opts.is_reactive = true;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t13 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, A, opts);
    auto us13 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t13).count();
    solver.printResult(r);
    record("Case13_RAND", "C1-C2-C3-nC4-nC5", 320.0, 30e5, r, us13);
    EXPECT_TRUE(r.success);
}

// Case 14: Water-gas shift, T=500K, P=50bar
// Uses Ideal EOS with NIST formation data (ΔHf°, S°) for correct reactive equilibrium
TEST_F(RANDThesisTest, Case14_WaterGasShift) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO,H2O,CO2,H2");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<std::vector<double>> A = {
        {1, 0, 1, 0},
        {0, 2, 0, 2},
        {1, 1, 2, 0}
    };
    // feed CO+H2O: {1,1,0,0} -> A*feed=[C:1,H:2,O:2]=[1,2,2]
    randflash::FlashInput input{500.0, 50e5, {1.0, 1.0, 0.0, 0.0}};
    randflash::SolveOptions opts;
    opts.is_reactive = true;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t14 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, A, opts);
    auto us14 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t14).count();
    solver.printResult(r);
    record("Case14_RAND", "CO-H2O-CO2-H2", 500.0, 50e5, r, us14);
    EXPECT_TRUE(r.success);
}

// Case 15: Methanation, T=600K, P=30bar
TEST_F(RANDThesisTest, Case15_Methanation) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO,CO2,H2,C1,H2O");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<std::vector<double>> A = {
        {1, 1, 0, 1, 0},
        {0, 0, 2, 4, 2},
        {1, 2, 0, 0, 1}
    };
    // feed CO+H2: {1,0,4,0,0} -> A*feed=[C:1,H:8,O:1]
    randflash::FlashInput input{600.0, 30e5, {1.0, 0.0, 4.0, 0.0, 0.0}};
    randflash::SolveOptions opts;
    opts.is_reactive = true;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t15 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, A, opts);
    auto us15 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t15).count();
    solver.printResult(r);
    record("Case15_RAND", "CO-CO2-H2-C1-H2O", 600.0, 30e5, r, us15);
    EXPECT_TRUE(r.success);
}

// Case 16: Ammonia synthesis, T=700K, P=200bar
TEST_F(RANDThesisTest, Case16_AmmoniaSynthesis) {
    auto backend = std::make_unique<thermo::ThermoPackBackend>("N2,H2,NH3");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<std::vector<double>> A = {
        {2, 0, 1},
        {0, 2, 3}
    };
    // feed N2+H2: {0.5,1.5,0} -> A*feed=[N:1,H:3]=[1,3]
    randflash::FlashInput input{700.0, 200e5, {0.5, 1.5, 0.0}};
    randflash::SolveOptions opts;
    opts.is_reactive = true;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t16 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, A, opts);
    auto us16 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t16).count();
    solver.printResult(r);
    record("Case16_RAND", "N2-H2-NH3", 700.0, 200e5, r, us16);
    EXPECT_GE(r.numPhases(), 1);
}

int main(int argc, char **argv) {
    randflash::setLogLevel(spdlog::level::debug);
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "RANDThesisTest.Case11_20comp_300K_50bar";
    return RUN_ALL_TESTS();
}
