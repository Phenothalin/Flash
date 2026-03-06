#!/usr/bin/env julia

"""
Gibbs Energy Verification using Clapeyron.jl

This script calculates Gibbs energies for ethane cracking system using Clapeyron
to compare with ThermoPack/SRK results.

System: C2H6 ⇌ C2H4 + H2
Temperature: 800 K
Pressure: 75000 Pa
"""

using Clapeyron

println("="^70)
println("Gibbs Energy Verification with Clapeyron.jl")
println("="^70)

# System setup
components = ["ethane", "ethylene", "hydrogen"]
T = 800.0  # K
P = 75000.0  # Pa

println("\nSystem Configuration:")
println("  Components: ", join(components, ", "))
println("  Temperature: $T K")
println("  Pressure: $P Pa")
println("  EOS: SRK (Soave-Redlich-Kwong)")

# Create SRK model
model = SRK(components)

println("\n" * "="^70)
println("Element Conservation Check")
println("="^70)

# Element matrix for ethane cracking
# A[element, component] = number of atoms
#        C2H6  C2H4  H2
# C  [    2     2     0  ]
# H  [    6     4     2  ]

A = [2.0 2.0 0.0;
     6.0 4.0 2.0]

# Feed: pure ethane [1, 0, 0] → elements: C=2, H=6
feed_moles = [1.0, 0.0, 0.0]
element_moles = A * feed_moles

println("Feed composition: ", feed_moles)
println("Element moles from feed: C=$(element_moles[1]), H=$(element_moles[2])")

# Helper function to check element conservation
function check_elements(n, name)
    elem = A * n
    total = sum(n)
    x = n ./ total

    println("\n$name:")
    println("  Mole numbers: [$(n[1]), $(n[2]), $(n[3])]")
    println("  Total moles: $total")
    println("  Composition x: [$(x[1]), $(x[2]), $(x[3])]")
    println("  Element moles: C=$(elem[1]), H=$(elem[2])")
    println("  Element error: ΔC=$(abs(elem[1] - element_moles[1])), ΔH=$(abs(elem[2] - element_moles[2]))")

    return elem
end

# Solution 1: Uniform distribution (from ThermoPack test)
n1 = [0.505567, 0.494433, 0.494433]

# Solution 2: User's proposed solution
n2 = [0.85 * 1.08145, 0.075 * 1.08145, 0.075 * 1.08145]

println("\n" * "="^70)
println("Solution 1: Uniform Distribution")
println("="^70)
check_elements(n1, "Solution 1")

println("\n" * "="^70)
println("Solution 2: User's Proposed Solution")
println("="^70)
check_elements(n2, "Solution 2")

# Helper function to calculate Gibbs energy
function calc_gibbs(n, name)
    total = sum(n)
    x = n ./ total

    println("\n$name:")
    println("  Composition x: [$(x[1]), $(x[2]), $(x[3])]")

    # Calculate Gibbs free energy using Clapeyron
    # Note: Clapeyron's gibbs_free_energy returns G in J
    G = gibbs_free_energy(model, P, T, n)

    println("  Total Gibbs energy G = $G J")

    return G
end

println("\n" * "="^70)
println("Gibbs Energy Comparison")
println("="^70)

G1 = calc_gibbs(n1, "Solution 1 (uniform)")
G2 = calc_gibbs(n2, "Solution 2 (user proposed)")

println("\n" * "="^70)
println("Result")
println("="^70)
println("ΔG = G2 - G1 = $(G2 - G1) J")

if G2 < G1
    println("\n*** Solution 2 has LOWER Gibbs energy by $(G1 - G2) J ***")
    println("This confirms the user's claim that [0.85, 0.075, 0.075] is thermodynamically more stable.")
else
    println("\nSolution 1 (uniform) has lower Gibbs energy by $(G2 - G1) J")
    println("This matches ThermoPack/SRK results.")
end

