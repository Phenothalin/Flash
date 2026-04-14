#include "RR_vl.hpp"
#include "thermo_backend.hpp"
#include "rr_logger.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>

static const std::string SWEEP_PATH = "/Users/madao/Desktop/2606/thesis_data/case8_sweep.csv";

TEST(Case8SweepTest, RR_CO2_C3_Sweep) {
    {
        std::ifstream chk(SWEEP_PATH);
        if (!chk.good()) {
            std::ofstream hdr(SWEEP_PATH);
            hdr << "Method,T_K,P_bar,Beta\n";
        }
    }
    std::ofstream sw(SWEEP_PATH, std::ios::app);
    sw << std::fixed << std::setprecision(6);

    thermo::ThermoPackBackend backend("CO2,C3", "PR");
    std::vector<double> z = {0.5, 0.5};
    for (double T = 200.0; T <= 375.0; T += 2.0) {
        for (double P_bar = 1.0; P_bar <= 85.0; P_bar += 2.0) {
            try {
                PTFlash f(P_bar*1e5, T, z, backend);
                f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
                sw << "RR," << T << "," << P_bar << "," << f.getVaporFraction() << "\n";
            } catch (...) {
                sw << "RR," << T << "," << P_bar << ",-1\n";
            }
        }
    }
    EXPECT_TRUE(true);
}

int main(int argc, char **argv) {
    rrflash::setLogLevel(spdlog::level::warn);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
