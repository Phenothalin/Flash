# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Flash is a C++17 thermodynamic flash calculation engine for multi-phase equilibrium computations. It solves Gibbs free energy minimization problems for multi-component, multi-phase systems using Newton-Raphson iteration with line search.

## Build System

### Prerequisites

Set this environment variable before building:

```bash
export VCPKG_ROOT=/path/to/vcpkg
```

Note: ThermoPack is now embedded in the project under `external/thermopack/`, so `THERMOPACK_DIR` is no longer required.

### Build Commands

```bash
# Configure and build
cd /Users/madao/code/Flash
mkdir -p build && cd build
cmake ..
make

# Build specific target
make <target_name>

# Clean build
rm -rf build && mkdir build && cd build && cmake .. && make
```

### Available Targets

- `thermo` - Thermodynamic backend library
- `phase_stability` - Phase stability analysis library
- `randflash` - Main flash solver library
- `compressibility_factor_test` - EOS validation test
- `tp_thermo_check` - ThermoPack integration test
- `2phase_tests` - Two-phase flash GTest suite
- `2phase_complex_test` - Complex two-phase scenarios
- `3phase_rand_test` - Three-phase flash test

### Running Tests

```bash
# From build directory
./Thermo/tests/compressibility_factor_test
./Thermo/tests/tp_thermo_check
./RAND/tests/2phase_tests
./RAND/tests/2phase_complex_test
./RAND/tests/3phase_rand_test
```

## Architecture

### Layered Design

The codebase follows a strict layered architecture with clear dependency flow:

```
RAND (Flash Solver)
    ↓ depends on
PhaseStability (Stability Analysis)
    ↓ depends on
Thermo (Thermodynamic Backend)
    ↓ depends on
External Libraries (ThermoPack, Eigen3)
```

### Embedded ThermoPack

ThermoPack is embedded in `external/thermopack/` with the following structure:

```
external/thermopack/
├── CMakeLists.txt          # Creates imported thermopack target
├── include/
│   └── cppThermopack/      # C++ headers (11 files)
│       ├── cubic.h
│       ├── thermo.h
│       └── ...
└── lib/
    └── macos/              # Platform-specific libraries
        ├── libthermopack.dylib
        └── libthermopack.a
```

Build options:
- `THERMOPACK_USE_STATIC=OFF` (default): Link dynamically
- `THERMOPACK_USE_STATIC=ON`: Link statically

### Core Modules

**Thermo/** - Thermodynamic property calculations
- `thermo_backend.hpp`: Abstract `IThermoBackend` interface defining thermodynamic operations
  - `chemicalPotentials()`, `dmu_dn()`: Chemical potential and derivatives
  - `fugacityCoefficients()`, `lnFugacityCoefficients()`: Fugacity calculations
  - `vaporPhaseFlag()`, `liquidPhaseFlag()`: Phase type identifiers
  - `minGibbsPhaseFlag()`: Auto-select Gibbs-minimum root (ThermoPack `Phase::mingibbs`)
  - `compressibilityFactor()`: Z = PV/(nRT) for phase type determination
- `thermopack_adapter.hpp`: `ThermoAdapterTP` wraps ThermoPack's Cubic EOS (PR, SRK)
- Provides chemical potentials (via `chemical_potential_tv`, ideal + residual), fugacity coefficients, and their derivatives
- Chemical potential calculation: TP → TV conversion using `specific_volume`, then `chemical_potential_tv` with `PropertyFlag::total`

**PhaseStability/** - Phase stability analysis
- `phase_stability.hpp/cpp`: `PhaseStabilityAnalyzer` class
- Implements Tangent Plane Distance (TPD) minimization via Michelsen successive substitution
- Multi-seed approach to detect incipient phases
- Supports dual reference phases (vapor-like and liquid-like)

**RAND/** - Multi-phase flash solver
- `rand_flash.hpp/cpp`: Main `RandFlash` class with Newton-Raphson solver
- `rand_init.cpp`: Initialization routines for 2-phase, 3-phase, and N-phase systems
- `rand_solver.cpp`: High-level solver entry points with `SolveOptions` control
- `linear_solver.hpp`: Abstract interface for linear algebra operations
- `eigen.cpp`: Eigen-based implementation of linear solver

### Key Design Patterns

1. **Dependency Injection**: `RandFlash` receives `IThermoBackend` and `LinearSolverInterface` references
2. **Strategy Pattern**: Multiple initialization strategies based on phase count
3. **Interface Segregation**: Abstract interfaces decouple implementations

## Key Data Structures

### FlashInput
Defines the input for flash calculations:
- `z`: Overall composition (mole fractions)
- `T`: Temperature
- `P`: Pressure
- `phase_labels`: Phase identifiers (e.g., "vapor", "liquid")

### MultiFlashResult
Contains the solution:
- `phases`: Vector of `PhaseContext`, each containing:
  - `state`: `PhaseState` with T, P, moleNumbers, phaseFlag
  - `x`: Composition (mole fractions)
  - `mu`: Chemical potentials
  - `m`, `M`: Jacobian matrices
- `success`: Convergence status
- `iterations`: Number of iterations
- `mu_infinity_norm`: Convergence error metric
- `pressure`, `temperature`: System conditions
- Convenience methods: `beta(j)`, `n_phase(j)`, `numPhases()`

### SolveOptions
Controls solver behavior:
- `max_iterations`: Maximum Newton iterations
- `tolerance`: Convergence tolerance
- `enable_stability_analysis`: Toggle automatic phase detection
- `verbose`: Debug output control

## Algorithm Overview

### Phase Stability Analysis
1. Compute TPD for trial compositions using successive substitution
2. Test multiple random seeds to find all incipient phases
3. Use both vapor-like and liquid-like reference phases
4. Return list of stable phases with their compositions

### Multi-Phase Flash Solver
1. Initialize phase fractions and compositions (via `rand_init.cpp`)
2. Set all phases to use `minGibbsPhaseFlag` (ThermoPack auto-selects stable root)
3. Newton-Raphson iteration on reduced chemical potentials:
   - Compute Hessian and gradient
   - Fix Hessian positive-definiteness via tangent space projection
   - Solve linear system for Newton step
   - Apply line search with alpha stepping
4. Check convergence based on reduced chemical potential differences
5. Post-convergence: Determine actual phase types via compressibility factor (Z > 0.5 → Vapor)
6. Reorder phases: Vapor → Oil-like liquid → Water-rich liquid
7. Return phase fractions, compositions, and convergence status

### Reactive Systems Support

The RAND algorithm supports reactive equilibrium through element conservation:

**Element Matrix (A)**: Defines elemental composition of each species
- A[e][i] = number of atoms of element e in species i
- Dimensions: E (elements) × C (species)
- For non-reactive systems: A = identity matrix
- For reactive systems: A = actual chemical formula matrix

**Conservation Constraints**:
- Element conservation: ∑_j ∑_i A[e][i] × n_i^(j) = b_e (constant)
- Replaces species conservation in non-reactive systems
- Automatically satisfies reaction equilibrium via element potentials

**Example**: Hydrocarbon system (C1, C2, C3)
```
       C1   C2   C3
  C  [  1    2    3  ]   (carbon atoms)
  H  [  4    6    8  ]   (hydrogen atoms)
```

## Code Navigation Tips

### Entry Points
- `RandFlash::solve()` in `RAND/src/rand_solver.cpp` - Main multi-phase flash solver
- `PhaseStabilityAnalyzer::analyze()` in `PhaseStability/src/phase_stability.cpp` - Stability testing
- `ThermoAdapterTP` in `Thermo/include/thermopack_adapter.hpp` - Thermodynamic calculations

### Understanding the Solver
1. Start with `RAND/include/rand_flash.hpp` for public API
2. Read `RAND/tests/2phase_tests.cpp` for usage examples
3. Study `RAND/src/rand_flash.cpp` for Newton-Raphson core algorithm
4. Review `RAND/src/rand_init.cpp` for initialization strategies

### Adding New Features
- New thermodynamic models: Implement `IThermoBackend` interface
- New linear solvers: Implement `LinearSolverInterface` interface
- New initialization strategies: Add to `rand_init.cpp`
- New convergence criteria: Modify `checkConvergence()` in `rand_flash.cpp`

### Using Reactive Systems

**Element Matrix Builder** (`RAND/include/element_matrix_builder.hpp`):

```cpp
#include "element_matrix_builder.hpp"

// Define species by chemical formula
std::vector<SpeciesFormula> species = {
    parseFormula("C1", "CH4"),    // Methane
    parseFormula("C2", "C2H6"),   // Ethane
    parseFormula("C3", "C3H8")    // Propane
};

// Build element matrix
std::vector<std::string> elementNames;
auto A = buildElementMatrix(species, elementNames);
// Returns: A[0] = [1, 2, 3] (carbon), A[1] = [4, 6, 8] (hydrogen)
// elementNames = ["C", "H"]

// Compute element moles from feed composition
std::vector<double> z = {0.5, 0.3, 0.2};
auto elementMoles = computeElementMoles(A, z);
// Returns: [1.7, 5.4] (carbon and hydrogen moles)
```

**Running Flash with Element Matrix**:

```cpp
// Setup backend and solver
auto backend = std::make_unique<ThermoPackBackend>("C1,C2,C3", "PR");
auto linSolver = ls::createEigenSolver();
RandFlash flash(*backend, *linSolver);

// Prepare input
FlashInput input;
input.temperature = 300.0;  // K
input.pressure = 1e5;       // Pa
input.feedMoles = {0.5, 0.3, 0.2};

// IMPORTANT: For reactive systems, MUST disable stability test
// and manually specify phase count
SolveOptions options;
options.enable_stability_test = false;  // Stability analysis not applicable
options.forced_phase_count = 1;         // Manually specify number of phases
auto result = flash.solve(input, A, options);

// Element conservation is automatically enforced
// Verify: computeElementMoles(A, result.phases[j].state.moleNumbers) = elementMoles
```

**Important Notes for Reactive Systems**:
- Species moles are NOT conserved; only element moles are conserved
- Phase stability analysis is NOT applicable (species composition changes via reactions)
- Must use `enable_stability_test = false` and `forced_phase_count`
- Some species may not exist initially (generated through reactions)

## Current Development Status

**Branch**: `feature/reaction`

Recent changes:
- Embedded ThermoPack into project (`external/thermopack/`), eliminating external dependency
- Renamed `thermo_adapter.hpp` to `thermopack_adapter.hpp` for clarity
- Replaced assembled chemical potential with ThermoPack's native `chemical_potential_tv` (ideal + residual)
- Added switchable stability analysis control via `SolveOptions::enable_stability_analysis`
- Unified multi-phase interface (removed separate `solveMultiPhase` methods)
- Improved convergence criteria using reduced chemical potentials
- Separated initialization logic into `rand_init.cpp`
- Added `minGibbsPhaseFlag()` to auto-select Gibbs-minimum root during iteration
- Streamlined `MultiFlashResult` - removed redundant `n_phase`/`beta` fields, data now accessed via `phases`
- Post-convergence phase type determination using compressibility factor Z
- **NEW**: Extended RAND algorithm to support reactive systems via element conservation
- **NEW**: Added element matrix builder (`element_matrix_builder.hpp/cpp`) for chemical formula parsing
- **NEW**: Implemented reactive system initialization functions (`initializeReactiveSinglePhase`, `initializeReactiveMultiPhase`)
- **NEW**: Added comprehensive test suite for reactive systems (`reactive_flash_test.cpp`)
- **NEW**: Element conservation automatically enforced throughout Newton-Raphson iteration

## Dependencies

- **ThermoPack**: Thermodynamic property calculations (embedded in `external/thermopack/`)
- **Eigen3**: Linear algebra library (installed via vcpkg)
- **GTest**: Unit testing framework (installed via vcpkg)
- **vcpkg**: C++ package manager (must set `VCPKG_ROOT`)

## File Naming Conventions

- Headers: `.hpp` extension
- Implementation: `.cpp` extension
- Tests: `*_test.cpp` or `*_tests.cpp` suffix
- Each module has `include/`, `src/`, and `tests/` subdirectories
