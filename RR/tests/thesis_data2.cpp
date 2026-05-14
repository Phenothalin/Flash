#include "RR_vl.hpp"
#include "RR_vll.hpp"
#include "thermo_backend.hpp"
#include "rr_logger.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>

static const std::string CSV_PATH = "/Users/madao/Desktop/2606/thesis_data/rr_results.csv";
static const std::string LOG_PATH = "/Users/madao/Desktop/2606/thesis_data/rr_detailed.log";
static const std::string ITER_CSV = "/Users/madao/Desktop/2606/thesis_data/rr_iter_history.csv";
static const std::string INNER_ITER_CSV = "/Users/madao/Desktop/2606/thesis_data/rr_inner_beta_history.csv";

class RRThesisTest : public ::testing::Test {
protected:
    std::ofstream csv;

    void printInnerBetaHistory(const std::string& cas,
                               const std::vector<double>& inner_beta_history,
                               const std::vector<int>& inner_start_indices,
                               const std::vector<int>& inner_counts) {
        if (inner_start_indices.size() != inner_counts.size()) return;
        for (size_t outer = 0; outer < inner_counts.size(); ++outer) {
            int start = inner_start_indices[outer];
            int count = inner_counts[outer];
            for (int local = 0; local < count; ++local) {
                int global = start + local;
                if (global >= 0 &&
                    global < static_cast<int>(inner_beta_history.size())) {
                    std::cout << "INNER_BETA " << cas
                              << " outer=" << (outer + 1)
                              << " local=" << (local + 1)
                              << " global=" << (global + 1)
                              << " beta_res=" << inner_beta_history[global]
                              << "\n";
                }
            }
        }
    }

    void SetUp() override {
        std::ifstream check(CSV_PATH);
        bool exists = check.good();
        check.close();
        csv.open(CSV_PATH, std::ios::app);
        if (!exists)
            csv << "Case,System,T_K,P_bar,Method,Beta,Iterations,Time_us,Status\n";

        std::ifstream check2(ITER_CSV);
        bool iter_exists = check2.good();
        check2.close();
        std::ofstream iter_csv(ITER_CSV, std::ios::app);
        if (!iter_exists)
            iter_csv << "Case,Iteration,KError\n";

        std::ifstream check3(INNER_ITER_CSV);
        bool inner_exists = check3.good();
        check3.close();
        std::ofstream inner_csv(INNER_ITER_CSV, std::ios::app);
        if (!inner_exists)
            inner_csv << "Case,OuterIteration,InnerIterationGlobal,InnerIterationLocal,BetaResidual\n";
    }

    void TearDown() override { csv.close(); }

    void log2Phase(const std::string& cas, const std::string& sys, double T, double P,
                   const std::string& method, double beta, int iters, long time_us,
                   const std::vector<double>& y, const std::vector<double>& x,
                   const std::vector<double>& iter_history = {},
                   const std::vector<double>& inner_beta_history = {},
                   const std::vector<int>& inner_start_indices = {},
                   const std::vector<int>& inner_counts = {}) {
        csv << cas << "," << sys << "," << T << "," << P/1e5 << ","
            << method << "," << beta << "," << iters << "," << time_us << ",CONVERGED\n";
        // Write iteration history
        if (!iter_history.empty()) {
            std::ofstream iter_csv(ITER_CSV, std::ios::app);
            for (size_t i = 0; i < iter_history.size(); ++i)
                iter_csv << cas << "," << (i+1) << "," << iter_history[i] << "\n";
        }
        if (!inner_beta_history.empty() &&
            inner_start_indices.size() == inner_counts.size()) {
            std::ofstream inner_csv(INNER_ITER_CSV, std::ios::app);
            for (size_t outer = 0; outer < inner_counts.size(); ++outer) {
                int start = inner_start_indices[outer];
                int count = inner_counts[outer];
                for (int local = 0; local < count; ++local) {
                    int global = start + local;
                    if (global >= 0 &&
                        global < static_cast<int>(inner_beta_history.size())) {
                        inner_csv << cas << "," << (outer + 1) << "," << (global + 1)
                                  << "," << (local + 1) << ","
                                  << inner_beta_history[global] << "\n";
                    }
                }
            }
        }
        std::ofstream log(LOG_PATH, std::ios::app);
        log << "\n=== " << cas << " [" << method << "] ===\n"
            << "T=" << T << "K P=" << P/1e5 << "bar beta=" << beta << "\n"
            << "Vapor: ";
        for (size_t i = 0; i < y.size(); ++i) log << y[i] << (i+1<y.size()?", ":"");
        log << "\nLiquid: ";
        for (size_t i = 0; i < x.size(); ++i) log << x[i] << (i+1<x.size()?", ":"");
        log << "\n";
    }

    void log3Phase(const std::string& cas, const std::string& sys, double T, double P,
                   const rr_vll::VLLResult& r, long time_us = 0) {
        csv << cas << "," << sys << "," << T << "," << P/1e5 << ",VLL,"
            << r.vapor_fraction << ";" << r.liquid1_fraction << ";" << r.liquid2_fraction
            << "," << r.iterations << "," << time_us << "," << (r.converged?"CONVERGED":"FAILED") << "\n";
        std::ofstream log(LOG_PATH, std::ios::app);
        log << "\n=== " << cas << " ===\n"
            << "T=" << T << "K P=" << P/1e5 << "bar\n"
            << "Vapor(beta=" << r.vapor_fraction << "): ";
        for (size_t i = 0; i < r.vapor_comp.size(); ++i) log << r.vapor_comp[i] << (i+1<r.vapor_comp.size()?", ":"");
        log << "\nLiquid1(beta=" << r.liquid1_fraction << "): ";
        for (size_t i = 0; i < r.liquid1_comp.size(); ++i) log << r.liquid1_comp[i] << (i+1<r.liquid1_comp.size()?", ":"");
        log << "\nLiquid2(beta=" << r.liquid2_fraction << "): ";
        for (size_t i = 0; i < r.liquid2_comp.size(); ++i) log << r.liquid2_comp[i] << (i+1<r.liquid2_comp.size()?", ":"");
        log << "\n";
    }
};

// Case 1: CH4-C2H6 Newton vs Halley, T=230K, P=20bar
TEST_F(RRThesisTest, Case1_Newton_vs_Halley) {
    thermo::ThermoPackBackend backend("C1,C2", "SRK");
    double T = 230.0, P = 20e5;
    std::vector<double> z = {0.5, 0.5};

    PTFlash fn(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    fn.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    auto t1 = std::chrono::high_resolution_clock::now();
    long tn = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    fn.printResult();
    log2Phase("Case1", "C1-C2", T, P, "Newton", fn.getVaporFraction(), fn.getIterations(), tn, fn.getVapComp(), fn.getLiqComp());

    PTFlash fh(P, T, z, backend);
    t0 = std::chrono::high_resolution_clock::now();
    fh.calculate(ConvergenceMethod::HALLEY);
    t1 = std::chrono::high_resolution_clock::now();
    long th = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    fh.printResult();
    log2Phase("Case1", "C1-C2", T, P, "Halley", fh.getVaporFraction(), fh.getIterations(), th, fh.getVapComp(), fh.getLiqComp());

    EXPECT_GE(fn.getVaporFraction(), 0.0);
    EXPECT_GE(fh.getVaporFraction(), 0.0);
}

// Case 2: CO2-C3H8 stable form, T=270K, P=30bar
TEST_F(RRThesisTest, Case2_CO2_C3H8_StableForm) {
    thermo::ThermoPackBackend backend("CO2,C3", "SRK");
    double T = 270.0, P = 30e5;
    std::vector<double> z = {0.3, 0.7};

    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc2 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    log2Phase("Case2", "CO2-C3", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc2, f.getVapComp(), f.getLiqComp());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

// Case 3: N2-CH4 low-T high-P, T=220K, P=50bar
TEST_F(RRThesisTest, Case3_N2_CH4_LowTemp) {
    thermo::ThermoPackBackend backend("N2,C1", "SRK");
    double T = 220.0, P = 50e5;
    std::vector<double> z = {0.4, 0.6};

    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc3 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    log2Phase("Case3", "N2-C1", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc3, f.getVapComp(), f.getLiqComp());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

// Case 4: H2O-CH4-nC6 three-phase VLL, T=300K, P=100bar
TEST_F(RRThesisTest, Case4_H2O_CH4_nC6_VLL) {
    thermo::ThermoPackBackend backend("H2O,C1,nC6", "SRK");
    double T = 298.15, P = 50e5;
    std::vector<double> z = {0.3, 0.3, 0.4};

    rr_vll::VLLFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto r = f.calculate();
    long tc4 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult(r);
    log3Phase("Case4", "H2O-C1-nC6", T, P, r, tc4);
    EXPECT_TRUE(r.converged);
}

// Case 5: CO2-H2O-nC10 three-phase VLL, T=320K, P=150bar
TEST_F(RRThesisTest, Case5_CO2_H2O_nC10_VLL) {
    thermo::ThermoPackBackend backend("CO2,H2O,nC10", "SRK");
    double T = 320.0, P = 150e5;
    std::vector<double> z = {0.2, 0.3, 0.5};

    rr_vll::VLLFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto r = f.calculate();
    long tc5 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult(r);
    log3Phase("Case5", "CO2-H2O-nC10", T, P, r, tc5);
    EXPECT_TRUE(r.converged);
}

// Cases 6,7,9,10,11: additional RR tests for cross-comparison
TEST_F(RRThesisTest, Case6_C1_C2_250K_30bar) {
    thermo::ThermoPackBackend backend("C1,C2", "SRK");
    double T = 250.0, P = 30e5;
    std::vector<double> z = {0.5, 0.5};
    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc6 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    log2Phase("Case6_RR", "C1-C2", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc6, f.getVapComp(), f.getLiqComp());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

TEST_F(RRThesisTest, Case7_CO2_C3_310K_70bar) {
    thermo::ThermoPackBackend backend("CO2,C3", "SRK");
    double T = 310.0, P = 70e5;
    std::vector<double> z = {0.5, 0.5};
    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc7 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    log2Phase("Case7_RR", "CO2-C3", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc7, f.getVapComp(), f.getLiqComp());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

TEST_F(RRThesisTest, Case9_5comp_280K_50bar) {
    thermo::ThermoPackBackend backend("C1,C2,C3,nC4,nC5", "SRK");
    double T = 280.0, P = 50e5;
    std::vector<double> z = {0.85, 0.08, 0.04, 0.02, 0.01};
    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc9 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    printInnerBetaHistory("Case9_RR", f.getIterBetaResidualHistory(), f.getOuterInnerStartIndices(), f.getOuterInnerCounts());
    log2Phase("Case9_RR", "C1-C2-C3-nC4-nC5", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc9, f.getVapComp(), f.getLiqComp(), f.getIterKHistory(), f.getIterBetaResidualHistory(), f.getOuterInnerStartIndices(), f.getOuterInnerCounts());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

TEST_F(RRThesisTest, Case10_11comp_295K_20bar) {
    thermo::ThermoPackBackend backend("N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7", "SRK");
    double T = 295.0, P = 20e5;
    std::vector<double> z = {0.003,0.015,0.55,0.14,0.12,0.05,0.045,0.03,0.025,0.012,0.01};
    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc10 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    printInnerBetaHistory("Case10_RR", f.getIterBetaResidualHistory(), f.getOuterInnerStartIndices(), f.getOuterInnerCounts());
    log2Phase("Case10_RR", "11comp", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc10, f.getVapComp(), f.getLiqComp(), f.getIterKHistory(), f.getIterBetaResidualHistory(), f.getOuterInnerStartIndices(), f.getOuterInnerCounts());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

TEST_F(RRThesisTest, Case11_20comp_300K_50bar) {
    thermo::ThermoPackBackend backend("N2,CO2,C1,C2,C3,iC4,nC4,iC5,nC5,nC6,nC7,nC8,nC9,nC10,nC11,nC12,nC13,nC14,nC15,nC16", "SRK");
    double T = 300.0, P = 50e5;
    std::vector<double> z = {0.01,0.02,0.70,0.10,0.06,0.03,0.02,0.015,0.01,0.008,0.006,0.004,0.003,0.002,0.001,0.001,0.001,0.001,0.001,0.007};
    PTFlash f(P, T, z, backend);
    auto t0 = std::chrono::high_resolution_clock::now();
    f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    long tc11 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-t0).count();
    f.printResult();
    printInnerBetaHistory("Case11_RR", f.getIterBetaResidualHistory(), f.getOuterInnerStartIndices(), f.getOuterInnerCounts());
    log2Phase("Case11_RR", "20comp", T, P, "Newton", f.getVaporFraction(), f.getIterations(), tc11, f.getVapComp(), f.getLiqComp(), f.getIterKHistory(), f.getIterBetaResidualHistory(), f.getOuterInnerStartIndices(), f.getOuterInnerCounts());
    EXPECT_GE(f.getVaporFraction(), 0.0);
}

// Case 8: moved to RR/tests/case8_sweep_test.cpp

int main(int argc, char **argv) {
    rrflash::setLogLevel(spdlog::level::debug);
    ::testing::InitGoogleTest(&argc, argv);
    // ::testing::GTEST_FLAG(filter) = "RRThesisTest.Case9_5comp_280K_50bar:RRThesisTest.Case10_11comp_295K_20bar:RRThesisTest.Case11_20comp_300K_50bar";
    return RUN_ALL_TESTS();
}
