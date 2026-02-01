#ifndef ELEMENT_MATRIX_BUILDER_HPP
#define ELEMENT_MATRIX_BUILDER_HPP

#include <string>
#include <vector>
#include <map>

namespace randflash {

/**
 * @brief Represents a chemical species with its elemental composition
 *
 * Example: CH4 (methane) -> name="CH4", atoms={{"C", 1}, {"H", 4}}
 */
struct SpeciesFormula {
    std::string name;
    std::map<std::string, int> atoms;  // Element symbol -> atom count

    SpeciesFormula() = default;
    SpeciesFormula(const std::string& n) : name(n) {}
    SpeciesFormula(const std::string& n, const std::map<std::string, int>& a)
        : name(n), atoms(a) {}
};

/**
 * @brief Builds element matrix A from species formulas
 *
 * The element matrix A has dimensions E × C where:
 * - E = number of elements
 * - C = number of species (components)
 * - A[e][i] = number of atoms of element e in species i
 *
 * @param species Vector of species formulas
 * @param elementNames Output vector of element names (ordered as matrix rows)
 * @return Element matrix A[E][C]
 *
 * Example:
 *   species = [CH4, C2H6, C3H8]
 *   Returns:
 *        CH4  C2H6  C3H8
 *   C  [  1    2     3   ]
 *   H  [  4    6     8   ]
 *   elementNames = ["C", "H"]
 */
std::vector<std::vector<double>> buildElementMatrix(
    const std::vector<SpeciesFormula>& species,
    std::vector<std::string>& elementNames);

/**
 * @brief Parses a chemical formula string into SpeciesFormula
 *
 * Supported formats:
 * - Simple: "CH4" -> {C:1, H:4}
 * - With numbers: "C2H6" -> {C:2, H:6}
 * - Multi-digit: "C10H22" -> {C:10, H:22}
 * - Multiple elements: "H2O" -> {H:2, O:1}
 *
 * @param name Species name (for identification)
 * @param formula Chemical formula string
 * @return SpeciesFormula with parsed elemental composition
 *
 * @throws std::invalid_argument if formula is malformed
 */
SpeciesFormula parseFormula(const std::string& name,
                            const std::string& formula);

/**
 * @brief Computes element mole numbers from species composition
 *
 * Calculates b_e = ∑_i A[e][i] * z_i for each element e
 *
 * @param elementMatrix Element matrix A[E][C]
 * @param speciesMoles Species mole numbers or mole fractions [C]
 * @return Element mole numbers [E]
 */
std::vector<double> computeElementMoles(
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& speciesMoles);

/**
 * @brief Verifies element conservation
 *
 * Checks if ∑_i A[e][i] * n_i ≈ b_e for all elements
 *
 * @param elementMatrix Element matrix A[E][C]
 * @param speciesMoles Species mole numbers [C]
 * @param elementMoles Expected element mole numbers [E]
 * @param tolerance Relative tolerance for comparison
 * @return true if conservation is satisfied within tolerance
 */
bool verifyElementConservation(
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& speciesMoles,
    const std::vector<double>& elementMoles,
    double tolerance = 1e-10);

} // namespace randflash

#endif // ELEMENT_MATRIX_BUILDER_HPP
