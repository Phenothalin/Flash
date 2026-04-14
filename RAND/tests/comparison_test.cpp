// comparison_test.cpp - Cases 23-26: detailed composition output for comparison
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include "thermo_backend.hpp"

static const std::string LOG_PATH = "/Users/madao/Desktop/2606/thesis_data/comparison_results.log";

class ComparisonTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;
    std::ofstream log;

    void SetUp() override {
        linSolver = ls::createEigenSolver();
        log.open(LOG_PATH, std::ios::app);
    }
    void TearDown() override { log.close(); }

    std::vector<std::vector<double>> eye(size_t n) {
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) A[i][i] = 1.0;
        return A;
    }
};

// Case 1: CH4-C2H6, T=230K, P=20bar - RAND result for comparison with RR
TEST_F(ComparisonTest, Case1_C1_C2_RAND) {
    log << "\n=== Case 1 RAND: CH4-C2H6 T=230K P=20bar ===\n";
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{230.0, 20e5, {0.5, 0.5}};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t0 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(2), opts);
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    solver.printResult(r);
    if (r.success && r.numPhases() == 2) {
        auto n0=r.n_phase(0); double b0=r.beta(0);
        auto n1=r.n_phase(1); double b1=r.beta(1);
        double y0=n0[0]/b0, y1=n0[1]/b0, x0=n1[0]/b1, x1=n1[1]/b1;
        log << "beta_V=" << b0 << " iters=" << r.iterations << " mu_err=" << r.mu_infinity_norm << " time=" << us << "us\n";
        log << "y(C1)=" << y0 << " y(C2)=" << y1 << "\n";
        log << "x(C1)=" << x0 << " x(C2)=" << x1 << "\n";
    }
    EXPECT_TRUE(r.success);
}

// Case 2: CO2-C3H8, T=270K, P=30bar - RAND result
TEST_F(ComparisonTest, Case2_CO2_C3_RAND) {
    log << "\n=== Case 2 RAND: CO2-C3H8 T=270K P=30bar ===\n";
    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO2,C3", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{338.0, 65e5, {0.5, 0.5}};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t0 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(2), opts);
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    solver.printResult(r);
    log << "phases=" << r.numPhases() << " iters=" << r.iterations << " mu_err=" << r.mu_infinity_norm << " time=" << us << "us\n";
    if (r.success && r.numPhases() == 2) {
        auto n0=r.n_phase(0); double b0=r.beta(0);
        auto n1=r.n_phase(1); double b1=r.beta(1);
        log << "beta_V=" << b0 << " y(CO2)=" << n0[0]/b0 << " y(C3)=" << n0[1]/b0 << "\n";
        log << "beta_L=" << b1 << " x(CO2)=" << n1[0]/b1 << " x(C3)=" << n1[1]/b1 << "\n";
    } else if (r.success && r.numPhases() == 1) {
        auto n=r.n_phase(0); double b=r.beta(0);
        log << "single phase beta=" << b << " x(CO2)=" << n[0]/b << " x(C3)=" << n[1]/b << "\n";
    }
    EXPECT_TRUE(r.success);
}

// Case 3: N2-CH4, T=220K, P=50bar - RAND result
TEST_F(ComparisonTest, Case3_N2_C1_RAND) {
    log << "\n=== Case 3 RAND: N2-CH4 T=220K P=50bar ===\n";
    auto backend = std::make_unique<thermo::ThermoPackBackend>("N2,C1", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{220.0, 50e5, {0.4, 0.6}};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t0 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(2), opts);
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    solver.printResult(r);
    log << "phases=" << r.numPhases() << " iters=" << r.iterations << " mu_err=" << r.mu_infinity_norm << " time=" << us << "us\n";
    if (r.success) {
        for (size_t p = 0; p < r.numPhases(); ++p) {
            auto n=r.n_phase(p); double b=r.beta(p);
            log << "Phase" << p << "(beta=" << b << "): N2=" << n[0]/b << " C1=" << n[1]/b << "\n";
        }
    }
    EXPECT_TRUE(r.success);
}

// Case 23 uses RR - remove PTFlash, use only RAND for comparison
// Case 23: CH4-C2H6 binary VLE, T=250K, P=30bar - detailed compositions + K values
TEST_F(ComparisonTest, Case23_C1_C2_BinaryVLE) {
    log << "\n=== Case 23: CH4-C2H6 Binary VLE T=250K P=30bar ===\n";

    // RAND solver
    auto backend = std::make_unique<thermo::ThermoPackBackend>("C1,C2", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{250.0, 30e5, {0.5, 0.5}};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t23 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(2), opts);
    auto us23 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t23).count();
    solver.printResult(r);
    if (r.success && r.numPhases() == 2) {
        auto n0 = r.n_phase(0); double b0 = r.beta(0);
        auto n1 = r.n_phase(1); double b1 = r.beta(1);
        double y0=n0[0]/b0, y1=n0[1]/b0, x0=n1[0]/b1, x1=n1[1]/b1;
        log << "[RAND] beta_V=" << b0 << " beta_L=" << b1
            << " mu_err=" << r.mu_infinity_norm << " iters=" << r.iterations << " time=" << us23 << "us\n";
        log << "  y(C1)=" << y0 << " y(C2)=" << y1 << " sum_y=" << y0+y1 << "\n";
        log << "  x(C1)=" << x0 << " x(C2)=" << x1 << " sum_x=" << x0+x1 << "\n";
        log << "  K(C1)=" << y0/x0 << " K(C2)=" << y1/x1 << "\n";
    }
    EXPECT_TRUE(r.success);
}

// Case 24: 11-comp natural gas, T=280K, P=50bar - K values and compositions
TEST_F(ComparisonTest, Case24_11comp_KValues) {
    log << "\n=== Case 24: 11-comp Natural Gas T=280K P=50bar ===\n";
    std::string comps = "N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7";
    std::vector<std::string> names = {"N2","CO2","C1","C2","C3","iC4","nC4","iC5","nC5","nC6","nC7"};
    std::vector<double> z = {0.003,0.015,0.55,0.14,0.12,0.05,0.045,0.03,0.025,0.012,0.01};

    auto backend = std::make_unique<thermo::ThermoPackBackend>(comps, "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{280.0, 50e5, z};
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    auto t24 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(11), opts);
    auto us24 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t24).count();
    solver.printResult(r);
    if (r.success && r.numPhases() >= 2) {
        auto n0 = r.n_phase(0); double b0 = r.beta(0);
        auto n1 = r.n_phase(1); double b1 = r.beta(1);
        log << "[RAND] beta_V=" << b0 << " beta_L=" << b1
            << " mu_err=" << r.mu_infinity_norm << " iters=" << r.iterations << " time=" << us24 << "us\n";
        double sy=0, sx=0;
        for (size_t i = 0; i < names.size(); ++i) {
            double yi=n0[i]/b0, xi=n1[i]/b1;
            sy+=yi; sx+=xi;
            log << "  " << names[i] << ": y=" << yi << " x=" << xi << " K=" << yi/xi << "\n";
        }
        log << "  sum_y=" << sy << " sum_x=" << sx << "\n";
    }
    EXPECT_TRUE(r.success);
}

// Case 25: H2O-CH4-nC6 VLLE, T=298K, P=50bar
TEST_F(ComparisonTest, Case25_H2O_CH4_nC6_VLLE) {
    log << "\n=== Case 25: H2O-CH4-nC6 VLLE T=298K P=50bar ===\n";
    std::vector<std::string> names = {"H2O","C1","nC6"};

    auto backend = std::make_unique<thermo::ThermoPackBackend>("H2O,C1,nC6", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{298.0, 50e5, {0.3, 0.5, 0.2}};
    randflash::SolveOptions opts;
    opts.phase_determination_strategy = "stable";
    auto t25 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, eye(3), opts);
    auto us25 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t25).count();
    solver.printResult(r);
    log << "[RAND-Stable] phases=" << r.numPhases() << " iters=" << r.iterations
        << " mu_err=" << r.mu_infinity_norm << " time=" << us25 << "us\n";
    double sum_beta = 0;
    for (size_t p = 0; p < r.numPhases(); ++p) {
        auto n = r.n_phase(p); double b = r.beta(p); sum_beta += b;
        log << "  Phase" << p << "(beta=" << b << "): ";
        for (size_t i = 0; i < names.size(); ++i)
            log << names[i] << "=" << n[i]/b << (i+1<names.size()?" ":"\n");
    }
    log << "  sum_beta=" << sum_beta << "\n";
    EXPECT_TRUE(r.success);
}

// Case 26: Water-gas shift reactive, T=600K, P=10bar
TEST_F(ComparisonTest, Case26_WaterGasShift_600K) {
    log << "\n=== Case 26: Water-Gas Shift T=600K P=10bar ===\n";
    std::vector<std::string> names = {"CO","H2O","CO2","H2"};
    std::vector<std::vector<double>> A = {
        {1, 0, 1, 0},
        {0, 2, 0, 2},
        {1, 1, 2, 0}
    };
    // feed CO+H2: {1,0,0,1} -> A*feed=[C:1,H:2,O:1]
    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO,H2O,CO2,H2", "SRK");
    randflash::RandFlash solver(*backend, *linSolver);
    randflash::FlashInput input{600.0, 10e5, {1.0, 0.0, 0.0, 1.0}};
    randflash::SolveOptions opts;
    opts.is_reactive = true;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 1;
    auto t26 = std::chrono::high_resolution_clock::now();
    auto r = solver.solve(input, A, opts);
    auto us26 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t26).count();
    solver.printResult(r);
    log << "[RAND] success=" << r.success << " iters=" << r.iterations
        << " mu_err=" << r.mu_infinity_norm << " time=" << us26 << "us\n";
    if (r.success) {
        auto n = r.n_phase(0);
        double total = 0; for (auto v : n) total += v;
        log << "  Moles: ";
        for (size_t i = 0; i < names.size(); ++i)
            log << names[i] << "=" << n[i] << (i+1<names.size()?" ":"\n");
        log << "  Mole fracs: ";
        for (size_t i = 0; i < names.size(); ++i)
            log << names[i] << "=" << n[i]/total << (i+1<names.size()?" ":"\n");
        // Element conservation check
        std::vector<double> b_elem = {1.0, 2.0, 1.0};
        for (size_t e = 0; e < A.size(); ++e) {
            double sum = 0;
            for (size_t i = 0; i < n.size(); ++i) sum += A[e][i]*n[i];
            log << "  Element" << e << " balance: " << sum << " (target=" << b_elem[e] << ")\n";
        }
    }
    EXPECT_TRUE(r.success);
}

int main(int argc, char **argv) {
    randflash::setLogLevel(spdlog::level::debug);
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::GTEST_FLAG(filter) = "ComparisonTest.Case2_CO2_C3_RAND";
    return RUN_ALL_TESTS();
}
