# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

LFT ("Lattice Field Theory") is a C++20 Monte Carlo code for lattice spin/field models: Ising (1D/2D/3D), φ⁴, XY, Potts and, the main current focus, the **Edwards–Anderson (EA) spin glass** with parallel tempering, multiple replicas and Houdayer cluster moves (2D and 3D). Results are analysed in R Markdown / Python notebooks under `notebooks/`.

## Build & test

CMake fetches every dependency with FetchContent (spdlog, Catch2, Lyra, CLI11, glm, yaml-cpp). OpenMP is required. On macOS, libomp is expected from Homebrew. The CLion build dirs `cmake-build-{debug,release,relwithdebinfo}` (Ninja) already exist and are gitignored.

```bash
cmake -S . -B cmake-build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build-release --target ea_pt_mt      # one target
cmake --build cmake-build-release --target all_targets   # all drivers, EA drivers and the test binary
```

Executables: `ising`, `ising1d`, `ising3d`, `phi4`, `xy` (from `src/Drivers/`); `ea_pt_mt`, `ea_pt_mt_3D` (from `src/EdwardsAnderson/`); `test`.

Tests use Catch2 and are all linked into a single `test` binary (`tests/`). They are not registered with CTest. When you add a test source, also add it to the `test` executable in `tests/CMakeLists.txt`. That target has `src/` on its include path, so EA headers are included as `"EdwardsAnderson/..."`.

The tests cover mostly the core lattice, field, Ising and sweep code. For EA, the only C++ test is `test_houdayer.cpp`. EA physics is otherwise checked against exact enumeration and other references in `notebooks/EdwardsAnderson/MC tests/`.

```bash
cmake --build cmake-build-debug --target test
./cmake-build-debug/tests/test                 # all tests
./cmake-build-debug/tests/test "ising sweep"   # one test case by name
./cmake-build-debug/tests/test "[Lattice]"     # by tag
./cmake-build-debug/tests/test --list-tests
```

`-O3 -ffast-math` is applied globally, so floating-point results are not bit-exact across changes.

## Architecture

### Core library: `src/LFT` (header-only INTERFACE target `lft`)
Include paths are relative to `src/LFT`, for example `#include "Field/Field.h"`. The executables also add `src/`, for example `#include "utils/fs.h"`.

- `Field/Lattice.h`: `MultiIndex<I,D>` and `Lattice<I=uint32_t, D=2>`. The lattice is periodic. It precomputes a flat nearest-neighbour table (`up(i,d)`, `dn(i,d)`) and the even/odd checkerboard site lists. Memory order is `'F'` (default) or `'C'`. The EA drivers use `'C'`.
- `Field/Field.h`: `Field<F, L>` is a flat `std::vector` over a lattice, with `corona` helpers. `write(ostream)` dumps the raw bytes, and that is the binary format of the `cfg_*.bin` files.
- `MonteCarlo/sweep.h`: generic `sweep` / `sweep_mt` that update even sites, then odd sites (checkerboard). The update is any callable `update(field, site[, rng])`.
- Model headers (`Ising/`, `Phi4/`, `XY/`, `Potts/`) define field types, heat-bath/Metropolis updaters and observables in per-model namespaces.

### Edwards–Anderson: `src/EdwardsAnderson` (namespace `lft::ea`)
- `ea.h`: `SpinField<L>` (int8 spins). The couplings live in `JField<F,L>`, defined on a lattice of dimension `DIM+1` whose leading extent is the direction. The coupling from site `i` in direction `d` is stored at `j[i + n_elements*d]`. This file also holds the `HeathBath` updater (log-odds acceptance), the J initialisers (`init_bernoulli`, `init_gaussian`) and the observables `energy`, `magnetisation`, `overlap` and `link_overlap`.
- `replicas.h`: `Replicas<L>` owns raw `SpinField*` pointers for the replicas at one β.
- `parallel_tempering_mt.h`: `ParallelTemperingMT<L>` holds `chains[beta_idx][replica_idx]`. Tempering swaps pointers between adjacent β (`exchange`). `sweep`, `exchange` and `houdayer` run in parallel over β chains with OpenMP. Each chain j has its own `taus` stream `rngs[j]` from a `utils/rand.h` `taus_array` with one stream per β, used for its sweeps and Houdayer moves and for the exchange of the pair (j, j+1), so results are identical for any `--n-threads`. Houdayer cluster moves are in `houdayer.h`.
- Driver: `ea_pt_mt.cpp` is one dimension-templated driver (`run<DIM>`). It is built twice, as `ea_pt_mt` (`EA_DIM=2`) and `ea_pt_mt_3D` (`EA_DIM=3`), so a change applies to both. Shared helpers are in `ea_io.h` (`init_j_field`, `init_spin_field`, `measure_em`), and `make_file_path(dir, prefix, dims, name, ext)` in `utils/fs.h` builds the file names for either dimension. One Monte Carlo step (sweep, exchange, Houdayer) is `mc_step`. Both drivers use CLI11 (`options_cli11.{h,cpp}`). The old single-thread `ea` driver and its Lyra options were removed (lft-d8t).

### EA run configuration and outputs
The CLI11 options accept a TOML config through `--config file.toml`. Keys use dashes, for example `j-file`, `n-replicas` and `beta = [...]`. Examples are in `EA_runs/*.toml`. CLI flags override the config.

Each run writes the following files into `--data-dir`. File names come from `utils/fs.cpp` `make_file_path`: `{prefix}_{Lx:03}x{Ly:03}_{name}.{ext}` in 2D and `{prefix}_{Lx:02}x{Ly:02}x{Lz:02}_{name}.{ext}` in 3D.
- `opt_*.yaml` / `opt_*.toml`: the effective options. The TOML file can be fed back with `--config`.
- `j_*.txt`: the couplings that were used, one value per line with 9 significant digits, so they read back bit-exactly. Pass it to `--j-file` to reuse the same disorder realisation. The reader requires exactly DIM·N values and rejects a file for another lattice size or dimension.
- `em_*_bNNN.txt`, one per β index, written with 17 significant digits: per replica, energy and magnetisation (totals, not per site). When there are 2 or more replicas, the q overlap and link overlap of replicas 0 and 1 follow.
- `cfg_*_bNNN.bin`: raw int8 configurations for all replicas, every `save-freq` sweeps.
- `exch_*.txt`: parallel-tempering acceptance statistics.

The notebooks read these files. See `notebooks/EdwardsAnderson/MC tests/test_utils.py` (`read_options`, `read_em_files`). Changing a file name or column layout breaks them.

### Conventions
- **Memory order is part of behaviour.** `Lattice` defaults to Fortran (`'F'`) order, while the simulation drivers explicitly use `'C'`. The indexing tests come in pairs: a default-order file and a `_C` file. Keep both paths working when you change indexing, neighbours, serialization or coordinate conventions.
- **Use the precomputed lattice structure.** `Lattice` precomputes the periodic neighbours and the checkerboard split, and sweeps update all even sites before the odd ones, including the OpenMP variants. Model updates should use `lat.up`, `lat.dn` or `Field::corona` rather than reimplementing boundary handling.
- **Where code lives.** Model code is template-based and lives in headers under `src/LFT`. Executable orchestration and non-template utilities live in `.cpp` files.
- **RNG handling.** RNGs are passed by reference into update kernels. Never share one mutable engine across OpenMP workers. The EA code gives each β chain its own stream (indexed by chain, not by `omp_get_thread_num()`), which keeps results independent of the thread count; keep it that way.
- **Raw configuration output.** `Field::write` writes contiguous binary field values with no embedded dimensions or metadata. Keep the element type, lattice order, `make_file_path` naming and the companion option files compatible with the analysis notebooks.
- **CLI libraries.** The Ising, φ⁴ and XY drivers parse options with Lyra. The multithreaded EA drivers use CLI11. To add an EA option:
  - declare it in `Options`;
  - register it in `Options::Options()` (`options_cli11.cpp`);
  - add it to `Options::emit()` if it should be recorded with the run.

### Issue tracking
Known bugs and cleanup tasks from a full code review are tracked with beads (`bd`, prefix `lft-`). Use `bd ready` to list them, `bd show <id>` for details, and `bd close <id> --reason ...` to close one. The issue database syncs to GitHub (`refs/dolt/data`) with `bd dolt push` / `bd dolt pull`.

### Notebooks
`jupytext.toml` pairs `ipynb` with `Rmd`. Only the `.Rmd` files are tracked, and `*.ipynb` is gitignored, so edit the `.Rmd` files. `scripts/run_*_jobs.py` and `ising/run.py` are Python drivers that call the binaries in `cmake-build-release/`.

## Open items (from TODO.adoc)
- Clean up the code and documentation.
- Parallelize EA at the level of the parallel-tempering chains.
- 3D EA model.
