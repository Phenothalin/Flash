#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include "thermo_backend.hpp"

static const std::string SWEEP_PATH = "/Users/madao/Desktop/2606/thesis_data/case8_sweep.csv";

class Case8SweepTest : public ::testing::Test {
protected:
    std::shared_ptr<ls::LinearSolverInterface> linSolver;
    void SetUp() override { linSolver = ls::createEigenSolver(); }

    std::vector<std::vector<double>> eye(size_t n) {
        std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
        for (size_t i = 0; i < n; ++i) A[i][i] = 1.0;
        return A;
    }
};

TEST_F(Case8SweepTest, RAND_CO2_C3_Sweep) {
    {
        std::ifstream chk(SWEEP_PATH);
        if (!chk.good()) {
            std::ofstream hdr(SWEEP_PATH);
            hdr << "Method,T_K,P_bar,Phases,Beta\n";
        }
    }
    std::ofstream sw(SWEEP_PATH, std::ios::app);
    sw << std::fixed << std::setprecision(6);

    auto backend = std::make_unique<thermo::ThermoPackBackend>("CO2,C3", "PR");
    randflash::RandFlash solver(*backend, *linSolver);
    std::vector<double> z = {0.5, 0.5};
    auto A = eye(2);
    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};
    for (double T = 200.0; T <= 375.0; T += 2.0) {
        for (double P_bar = 1.0; P_bar <= 85.0; P_bar += 2.0) {
            randflash::FlashInput input{T, P_bar*1e5, z};
            auto r = solver.solve(input, A, opts);
            double beta_v = -1.0;
            if (r.success) {
                beta_v = 0.0;
                for (size_t p = 0; p < r.numPhases(); ++p)
                    if (r.phases[p].state.phaseFlag == backend->vaporPhaseFlag())
                        beta_v += r.beta(p);
            }
            sw << "RAND," << T << "," << P_bar << "," << r.numPhases() << "," << beta_v << "\n";
        }
    }
    EXPECT_TRUE(true);
}

int main(int argc, char **argv) {
    randflash::setLogLevel(spdlog::level::warn);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
