#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include "rand_flash.hpp"
#include "rand_logger.hpp"
#include "thermo_backend.hpp"

static const std::string TXY_PATH = "/Users/madao/Desktop/2606/thesis_data/txy_c2c3_rand.csv";

TEST(TxyC2C3Test, RAND_C2C3_Txy) {
    std::ofstream out(TXY_PATH);
    out << "z_C2,T_bubble_K,T_dew_K,y_C2_bubble,x_C2_dew\n";
    out << std::fixed << std::setprecision(6);

    auto backend = std::make_unique<thermo::ThermoPackBackend>("C2,C3", "SRK");
    auto linSolver = ls::createEigenSolver();
    randflash::RandFlash solver(*backend, *linSolver);

    const double P = 10.0e5; // 10 bar
    std::vector<std::vector<double>> A = {{1,0},{0,1}};

    randflash::SolveOptions opts;
    opts.enable_stability_test = false;
    opts.forced_phase_count = 2;
    opts.forced_phase_flags = {backend->vaporPhaseFlag(), backend->liquidPhaseFlag()};

    for (int i = 0; i <= 50; ++i) {
        double z_C2 = i * 0.02;
        std::vector<double> z = {z_C2, 1.0 - z_C2};

        // Bubble point: feed = liquid, sweep T from 240 to 310 K
        double T_bubble = -1.0, y_C2 = -1.0;
        for (double T = 240.0; T <= 310.0; T += 0.5) {
            randflash::FlashInput input{T, P, z};
            auto r = solver.solve(input, A, opts);
            if (r.success) {
                double beta_v = 0.0;
                int vi = -1;
                for (size_t p = 0; p < r.numPhases(); ++p)
                    if (r.phases[p].state.phaseFlag == backend->vaporPhaseFlag()) {
                        beta_v = r.beta(p); vi = (int)p;
                    }
                if (beta_v > 0.01 && beta_v < 0.99) {
                    T_bubble = T;
                    if (vi >= 0) y_C2 = r.phases[vi].x[0];
                    break;
                }
            }
        }

        // Dew point: feed = vapor, sweep T from 310 to 240 K
        double T_dew = -1.0, x_C2 = -1.0;
        for (double T = 310.0; T >= 240.0; T -= 0.5) {
            randflash::FlashInput input{T, P, z};
            auto r = solver.solve(input, A, opts);
            if (r.success) {
                double beta_v = 0.0;
                int li = -1;
                for (size_t p = 0; p < r.numPhases(); ++p) {
                    if (r.phases[p].state.phaseFlag == backend->vaporPhaseFlag())
                        beta_v = r.beta(p);
                    else
                        li = (int)p;
                }
                if (beta_v > 0.01 && beta_v < 0.99) {
                    T_dew = T;
                    if (li >= 0) x_C2 = r.phases[li].x[0];
                    break;
                }
            }
        }

        out << z_C2 << "," << T_bubble << "," << T_dew << "," << y_C2 << "," << x_C2 << "\n";
    }
    EXPECT_TRUE(true);
}

int main(int argc, char **argv) {
    randflash::setLogLevel(spdlog::level::warn);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
