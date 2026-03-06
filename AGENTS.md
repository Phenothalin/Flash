# Repository Guidelines

## Project Structure & Module Organization
Flash is a C++17 multi-module thermodynamic flash engine.
- `RAND/`: main multi-phase solver (`include/`, `src/`, `tests/`)
- `RR/`: Rachford-Rice flash solver (`include/`, `src/`, `tests/`)
- `PhaseStability/`: Michelsen TPD phase-stability analysis
- `Thermo/`: ThermoPack backend adapter and thermo checks
- `external/thermopack/`: embedded ThermoPack headers and prebuilt libs
- `build/`: out-of-source build artifacts (do not commit generated files)

Keep dependencies layered: `RAND -> PhaseStability -> Thermo -> external`.

## Build, Test, and Development Commands
Prerequisite: set `VCPKG_ROOT` and use CMake 3.18+.

```bash
# Configure (from repo root)
cmake -S . -B build -G Ninja

# Build all targets
cmake --build build

# Build one target (example)
cmake --build build --target randflash
```

Run tests from `build/` by executable target name, for example:
```bash
./RAND/2phase_tests
./RAND/stable_method_test
./RR/RRflash_test
./Thermo/compressibility_factor_test
```
On Windows, use corresponding `.exe` paths (see `GUIDE.md`).

## Coding Style & Naming Conventions
- Language standard: C++17 (`CMAKE_CXX_STANDARD 17`).
- Indentation: 2 spaces; braces on same line (follow existing files).
- Types/classes: `PascalCase` (for example `RandFlash`, `SolveOptions`).
- Functions/variables/files: `snake_case` (for example `solve_reactive`, `rand_solver.cpp`).
- Headers in `include/`, implementations in `src/`; keep module APIs minimal and explicit.
- No enforced formatter config is checked in; match nearby code style before submitting.

## Testing Guidelines
- Framework: GoogleTest (`GTest::gtest`, `GTest::gtest_main`).
- Add tests under each module’s `tests/` directory.
- Prefer deterministic fixture-based tests with `TEST_F` for solver scenarios.
- Name test files by feature, such as `reactive_flash_test.cpp` or `RR_vll_test.cpp`.
- Run affected module tests locally before opening a PR.

## Commit & Pull Request Guidelines
- Commit messages in history are short, action-oriented, and often scoped by module.
- Use concise subject lines; optional prefixes like `fix:`, `chore:`, `docs:` are acceptable.
- For multi-part changes, use numbered bullets in the body.
- PRs should include: purpose, touched modules, test evidence (commands + key output), and platform notes (macOS/Windows) when relevant.
