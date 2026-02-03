#include "rand_utils.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace randflash {
namespace utils {

std::vector<double> normalize(
    const std::vector<double>& v,
    double min_floor)
{
    std::vector<double> x = v;

    // Apply minimum floor to all elements
    for (double& xi : x) {
        xi = std::max(xi, min_floor);
    }

    // Calculate sum
    double s = std::accumulate(x.begin(), x.end(), 0.0);

    // If sum is invalid, return uniform distribution
    if (s <= 0.0) {
        const double uni = 1.0 / static_cast<double>(x.size());
        std::fill(x.begin(), x.end(), uni);
    } else {
        // Normalize to sum = 1.0
        for (double& xi : x) {
            xi /= s;
        }
    }

    return x;
}

std::vector<double> buildWilsonGuess(
    const std::vector<double>& feedMoles,
    const std::vector<double>& K,
    bool vapor_like)
{
    std::vector<double> x = feedMoles;

    // Normalize feed first
    const double s = std::accumulate(x.begin(), x.end(), 0.0);
    if (s > 0.0) {
        for (double& xi : x) {
            xi /= s;
        }
    }

    // Apply K-values
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] = vapor_like ? (x[i] * K[i]) : (x[i] / K[i]);
    }

    // Normalize result
    return normalize(x);
}

int detectWaterComponent(
    const std::vector<std::string>& names)
{
    for (size_t i = 0; i < names.size(); ++i) {
        const std::string& name = names[i];
        if (name == "H2O" || name == "WATER" ||
            name == "Water" || name == "water") {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool isWaterSystem(
    const std::vector<std::string>& names,
    const std::vector<double>& feedMoles,
    double threshold)
{
    int water_idx = detectWaterComponent(names);
    if (water_idx < 0) {
        return false;
    }

    // Calculate water mole fraction
    double total = std::accumulate(feedMoles.begin(), feedMoles.end(), 0.0);
    if (total <= 0.0) {
        return false;
    }

    double water_frac = feedMoles[water_idx] / total;
    return water_frac > threshold;
}

} // namespace utils
} // namespace randflash
