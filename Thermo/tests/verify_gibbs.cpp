#include <iostream>
#include <vector>
#include <numeric>
#include "thermo_backend.hpp"
#include "element_matrix_builder.hpp"
#include "rand_flash.hpp"  // For PhaseState

using namespace randflash;
using namespace thermo;

int main() {
    // Setup backend
    ThermoPackBackend backend("C2,C2_1,H2", "SRK");

    double T = 800.0;  // K
    double P = 75000;  // Pa

    // Element matrix for ethane cracking
    std::vector<SpeciesFormula> species = {
        parseFormula("C2", "C2H6"),
        parseFormula("C2_1", "C2H4"),
        parseFormula("H2", "H2")
    };
    std::vector<std::string> elementNames;
    auto A = buildElementMatrix(species, elementNames);

    // Original feed: [1, 0, 0] -> element moles: C=2, H=6
    std::vector<double> feedMoles = {1.0, 0.0, 0.0};
    auto elementMoles = computeElementMoles(A, feedMoles);

    std::cout << "Element moles from feed: C=" << elementMoles[0]
              << ", H=" << elementMoles[1] << "\n\n";

    // Solution 1: Uniform distribution (current result)
    std::vector<double> n1 = {0.505567, 0.494433, 0.494433};  // From test output
    double total1 = std::accumulate(n1.begin(), n1.end(), 0.0);

    // Solution 2: User's proposed better solution
    std::vector<double> n2 = {0.85 * 1.08145, 0.075 * 1.08145, 0.075 * 1.08145};
    double total2 = std::accumulate(n2.begin(), n2.end(), 0.0);

    // Verify element conservation for both
    auto check_elements = [&](const std::vector<double>& n, const std::string& name) {
        std::vector<double> elem(2, 0.0);
        for (size_t e = 0; e < 2; ++e) {
            for (size_t i = 0; i < 3; ++i) {
                elem[e] += A[e][i] * n[i];
            }
        }
        double total = std::accumulate(n.begin(), n.end(), 0.0);
        std::cout << name << ":\n";
        std::cout << "  Mole numbers: [" << n[0] << ", " << n[1] << ", " << n[2] << "]\n";
        std::cout << "  Total moles: " << total << "\n";
        std::cout << "  Element moles: C=" << elem[0] << ", H=" << elem[1] << "\n";
        std::cout << "  Element error: ΔC=" << std::abs(elem[0] - elementMoles[0])
                  << ", ΔH=" << std::abs(elem[1] - elementMoles[1]) << "\n";
        return elem;
    };

    std::cout << "=== Solution 1: Uniform distribution ===" << std::endl;
    check_elements(n1, "Solution 1");

    std::cout << "\n=== Solution 2: User's proposed solution ===" << std::endl;
    check_elements(n2, "Solution 2");

    // Calculate Gibbs energy for both
    std::cout << "\n=== Gibbs Energy Comparison ===" << std::endl;

    auto calc_gibbs = [&](const std::vector<double>& n, const std::string& name) {
        std::vector<double> x(3);
        double total = std::accumulate(n.begin(), n.end(), 0.0);
        for (size_t i = 0; i < 3; ++i) x[i] = n[i] / total;

        // Create PhaseState
        PhaseState st;
        st.Temperature = T;
        st.Pressure = P;
        st.moleNumbers = n;
        st.phaseFlag = backend.minGibbsPhaseFlag();

        auto mu = backend.chemicalPotentials(st);

        // G = Σ n_i μ_i
        double G = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            G += n[i] * mu[i];
        }

        std::cout << name << ":\n";
        std::cout << "  Composition x: [" << x[0] << ", " << x[1] << ", " << x[2] << "]\n";
        std::cout << "  Chemical potentials μ: [" << mu[0] << ", " << mu[1] << ", " << mu[2] << "] J/mol\n";
        std::cout << "  Total Gibbs energy G = Σ n_i μ_i = " << G << " J\n";

        return G;
    };

    double G1 = calc_gibbs(n1, "Solution 1 (uniform)");
    double G2 = calc_gibbs(n2, "Solution 2 (user proposed)");

    std::cout << "\n=== Result ===" << std::endl;
    std::cout << "ΔG = G2 - G1 = " << (G2 - G1) << " J" << std::endl;

    if (G2 < G1) {
        std::cout << "\n*** USER IS CORRECT! ***" << std::endl;
        std::cout << "Solution 2 has LOWER Gibbs energy by " << (G1 - G2) << " J" << std::endl;
        std::cout << "This proves the algorithm is stuck in a LOCAL minimum, not the global minimum." << std::endl;
        std::cout << "This is a pure optimization problem, NOT related to missing ΔG_f°." << std::endl;
    } else {
        std::cout << "\nSolution 1 (uniform) has lower Gibbs energy." << std::endl;
    }

    return 0;
}
