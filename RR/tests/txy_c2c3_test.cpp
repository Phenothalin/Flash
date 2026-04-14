#include "RR_vl.hpp"
#include "thermo_backend.hpp"
#include "rr_logger.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>

static const std::string TXY_PATH = "/Users/madao/Desktop/2606/thesis_data/txy_c2c3_rr.csv";

TEST(TxyC2C3Test, RR_C2C3_Txy) {
    std::ofstream out(TXY_PATH);
    out << "z_C2,T_bubble_K,T_dew_K,y_C2_bubble,x_C2_dew\n";
    out << std::fixed << std::setprecision(6);

    thermo::ThermoPackBackend backend("C2,C3", "SRK");
    const double P = 10.0e5; // 10 bar

    for (int i = 0; i <= 50; ++i) {
        double z_C2 = i * 0.02;
        std::vector<double> z = {z_C2, 1.0 - z_C2};

        // Bubble point: first T where two-phase appears from below
        double T_bubble = -1.0, y_C2 = -1.0;
        for (double T = 240.0; T <= 310.0; T += 0.5) {
            try {
                PTFlash f(P, T, z, backend);
                f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
                double bv = f.getVaporFraction();
                if (bv > 0.01 && bv < 0.99) {
                    T_bubble = T;
                    y_C2 = f.getVapComp()[0];
                    break;
                }
            } catch (...) {}
        }

        // Dew point: first T where two-phase appears from above
        double T_dew = -1.0, x_C2 = -1.0;
        for (double T = 310.0; T >= 240.0; T -= 0.5) {
            try {
                PTFlash f(P, T, z, backend);
                f.calculate(ConvergenceMethod::NEWTON_RAPHSON);
                double bv = f.getVaporFraction();
                if (bv > 0.01 && bv < 0.99) {
                    T_dew = T;
                    x_C2 = f.getLiqComp()[0];
                    break;
                }
            } catch (...) {}
        }

        out << z_C2 << "," << T_bubble << "," << T_dew << "," << y_C2 << "," << x_C2 << "\n";
    }
    EXPECT_TRUE(true);
}

int main(int argc, char **argv) {
    rrflash::setLogLevel(spdlog::level::warn);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
