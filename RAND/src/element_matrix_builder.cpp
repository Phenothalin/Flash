#include "element_matrix_builder.hpp"
#include <algorithm>
#include <cctype>
#include <set>
#include <stdexcept>
#include <cmath>

namespace randflash {

SpeciesFormula parseFormula(const std::string& name,
                            const std::string& formula) {
    SpeciesFormula result(name);

    if (formula.empty()) {
        throw std::invalid_argument("Empty formula string");
    }

    size_t i = 0;
    while (i < formula.size()) {
        // Skip whitespace
        if (std::isspace(formula[i])) {
            ++i;
            continue;
        }

        // Parse element symbol (uppercase letter followed by optional lowercase)
        if (!std::isupper(formula[i])) {
            throw std::invalid_argument(
                "Expected uppercase letter at position " + std::to_string(i) +
                " in formula: " + formula);
        }

        std::string element;
        element += formula[i++];

        // Add lowercase letters to element symbol
        while (i < formula.size() && std::islower(formula[i])) {
            element += formula[i++];
        }

        // Parse count (optional, defaults to 1)
        int count = 0;
        while (i < formula.size() && std::isdigit(formula[i])) {
            count = count * 10 + (formula[i] - '0');
            ++i;
        }

        if (count == 0) {
            count = 1;  // Default count
        }

        // Add to atoms map
        result.atoms[element] += count;
    }

    if (result.atoms.empty()) {
        throw std::invalid_argument("No elements parsed from formula: " + formula);
    }

    return result;
}

std::vector<std::vector<double>> buildElementMatrix(
    const std::vector<SpeciesFormula>& species,
    std::vector<std::string>& elementNames) {

    if (species.empty()) {
        throw std::invalid_argument("Empty species list");
    }

    // Collect all unique elements
    std::set<std::string> elementSet;
    for (const auto& sp : species) {
        for (const auto& [element, count] : sp.atoms) {
            elementSet.insert(element);
        }
    }

    // Convert to sorted vector for consistent ordering
    elementNames.clear();
    elementNames.assign(elementSet.begin(), elementSet.end());
    std::sort(elementNames.begin(), elementNames.end());

    const size_t E = elementNames.size();  // Number of elements
    const size_t C = species.size();        // Number of species

    // Build element matrix A[E][C]
    std::vector<std::vector<double>> A(E, std::vector<double>(C, 0.0));

    for (size_t i = 0; i < C; ++i) {
        for (size_t e = 0; e < E; ++e) {
            const std::string& element = elementNames[e];
            auto it = species[i].atoms.find(element);
            if (it != species[i].atoms.end()) {
                A[e][i] = static_cast<double>(it->second);
            }
        }
    }

    return A;
}

std::vector<double> computeElementMoles(
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& speciesMoles) {

    if (elementMatrix.empty()) {
        throw std::invalid_argument("Empty element matrix");
    }

    const size_t E = elementMatrix.size();       // Number of elements
    const size_t C = elementMatrix[0].size();    // Number of species

    if (speciesMoles.size() != C) {
        throw std::invalid_argument(
            "Species moles size (" + std::to_string(speciesMoles.size()) +
            ") does not match matrix columns (" + std::to_string(C) + ")");
    }

    std::vector<double> elementMoles(E, 0.0);

    for (size_t e = 0; e < E; ++e) {
        for (size_t i = 0; i < C; ++i) {
            elementMoles[e] += elementMatrix[e][i] * speciesMoles[i];
        }
    }

    return elementMoles;
}

bool verifyElementConservation(
    const std::vector<std::vector<double>>& elementMatrix,
    const std::vector<double>& speciesMoles,
    const std::vector<double>& elementMoles,
    double tolerance) {

    std::vector<double> computed = computeElementMoles(elementMatrix, speciesMoles);

    if (computed.size() != elementMoles.size()) {
        return false;
    }

    for (size_t e = 0; e < computed.size(); ++e) {
        double diff = std::abs(computed[e] - elementMoles[e]);
        double scale = std::max(std::abs(elementMoles[e]), 1.0);
        if (diff > tolerance * scale) {
            return false;
        }
    }

    return true;
}

} // namespace randflash
