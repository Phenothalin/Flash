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
- `beta`: Phase fractions
- `x`: Phase compositions (matrix: phases × components)
- `converged`: Convergence status
- `iterations`: Number of iterations
- `mu_reduced`: Reduced chemical potentials

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
2. Newton-Raphson iteration on reduced chemical potentials:
   - Compute Hessian and gradient
   - Fix Hessian positive-definiteness via tangent space projection
   - Solve linear system for Newton step
   - Apply line search with alpha stepping
3. Check convergence based on reduced chemical potential differences
4. Return phase fractions, compositions, and convergence status

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

## Current Development Status

**Branch**: `feature/multiphase`

Recent changes:
- Embedded ThermoPack into project (`external/thermopack/`), eliminating external dependency
- Renamed `thermo_adapter.hpp` to `thermopack_adapter.hpp` for clarity
- Replaced assembled chemical potential with ThermoPack's native `chemical_potential_tv` (ideal + residual)
- Added switchable stability analysis control via `SolveOptions::enable_stability_analysis`
- Unified multi-phase interface (removed separate `solveMultiPhase` methods)
- Improved convergence criteria using reduced chemical potentials
- Separated initialization logic into `rand_init.cpp`

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
