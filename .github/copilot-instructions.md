# Copilot instructions

## Build and test

This is a C++20 CMake project. OpenMP is required; the remaining C++ dependencies are downloaded with `FetchContent` during configuration.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

Tests are one Catch2 executable and are not registered with CTest:

```bash
cmake --build build --target test -j
./build/tests/test
./build/tests/test "Lattice indexing"  # one test case
./build/tests/test "[ising]"           # tests selected by tag
./build/tests/test --list-tests
```

When adding a test source, also add it to the `test` executable in `tests/CMakeLists.txt`.

## Architecture

- `src/LFT` is the reusable, mostly header-only simulation layer exposed through the `lft` INTERFACE target. `Field/Lattice.h` precomputes periodic nearest neighbors and checkerboard (even/odd) site lists; `Field/Field.h` couples storage to a lattice. Model headers under `Ising`, `Phi4`, `Potts`, and `XY` provide fields, observables, and update kernels, while `MonteCarlo/sweep.h` applies updates in even-then-odd order.
- `src/Drivers` contains standalone simulation programs (`ising`, `ising1d`, `ising3d`, `phi4`, and `xy`). They parse CLI options, construct the templated lattice/model types, run thermalization and measurement loops, and write data consumed by the R Markdown notebooks and Python scripts.
- `src/EdwardsAnderson` is a separate spin-glass application layer. It builds bond variables as a `(DIM + 1)`-dimensional `JField`, with the leading extent selecting bond direction. The parallel-tempering implementation organizes replicas by beta, performs local sweeps, neighboring-beta exchanges, and optional Houdayer cluster moves. `ea_pt_mt` and `ea_pt_mt_3D` use CLI11, TOML config files, OpenMP, and per-thread RNG state.
- `src/utils` supplies shared file naming, logging, RNG, hardware, and progress helpers used mainly by the executable drivers. Simulation outputs and run configurations live outside the library API; `EA_runs` contains representative TOML inputs.

## Repository conventions

- Lattice memory order is part of behavior. `Lattice` defaults to Fortran (`'F'`) order, while simulation drivers explicitly use C (`'C'`) order. Indexing tests intentionally have paired default-order and `_C` files; preserve both paths when changing indexing, neighbors, serialization, or coordinate conventions.
- Periodic neighbors and checkerboard partitions are precomputed by `Lattice`. Sweep implementations update all even sites before odd sites, including OpenMP variants; model updates should use `lat.up`, `lat.dn`, or `Field::corona` rather than reimplementing boundary handling.
- Core model code is template-based and normally lives in headers under `src/LFT`; executable orchestration and non-template utilities live in `.cpp` files. Public library includes are relative to `src/LFT` (for example, `"Field/Field.h"` and `"MonteCarlo/sweep.h"`).
- Random generators are passed by reference into update kernels. Multithreaded paths index a per-thread RNG collection with `omp_get_thread_num()`; do not share one mutable engine across OpenMP workers.
- Raw configuration output is written as contiguous binary field values without embedded dimensions or metadata. Keep element type, lattice order, output naming through `make_file_path`, and companion option files compatible with the analysis notebooks.
- The older Ising/Phi4/XY drivers use Lyra for CLI parsing. The newer multithreaded Edwards-Anderson drivers use CLI11 and support `--config` TOML files; add Edwards-Anderson options in `Options`, register them in `Options::Options()`, and include them in `Options::emit()` when they must be recorded with a run.
