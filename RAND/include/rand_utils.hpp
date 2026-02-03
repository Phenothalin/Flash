#pragma once

#include <vector>
#include <string>

namespace randflash {
namespace utils {

/**
 * @brief Normalize a composition vector
 *
 * Ensures all elements are >= min_floor and sum to 1.0.
 * If sum is zero or negative, returns uniform distribution.
 *
 * @param v Input vector to normalize
 * @param min_floor Minimum value for each element (default: 1e-14)
 * @return Normalized vector
 */
std::vector<double> normalize(
    const std::vector<double>& v,
    double min_floor = 1e-14);

/**
 * @brief Build Wilson K-value based composition guess
 *
 * Constructs a vapor-like or liquid-like composition from feed
 * using Wilson K-values: x_vapor = x_feed * K, x_liquid = x_feed / K
 *
 * @param feedMoles Feed composition (mole fractions or moles)
 * @param K Wilson K-values for each component
 * @param vapor_like If true, multiply by K; if false, divide by K
 * @return Normalized composition guess
 */
std::vector<double> buildWilsonGuess(
    const std::vector<double>& feedMoles,
    const std::vector<double>& K,
    bool vapor_like);

/**
 * @brief Detect water component index in component names
 *
 * Searches for common water identifiers: "H2O", "WATER", "Water", "water"
 *
 * @param names Vector of component names
 * @return Index of water component, or -1 if not found
 */
int detectWaterComponent(
    const std::vector<std::string>& names);

/**
 * @brief Check if system is water-dominated
 *
 * Returns true if water component exists and its mole fraction
 * exceeds the threshold.
 *
 * @param names Vector of component names
 * @param feedMoles Feed composition
 * @param threshold Minimum water mole fraction (default: 1e-4)
 * @return True if water-dominated system
 */
bool isWaterSystem(
    const std::vector<std::string>& names,
    const std::vector<double>& feedMoles,
    double threshold = 1e-4);

} // namespace utils
} // namespace randflash
