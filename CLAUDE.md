# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`educational_VMEC` is a heavily stripped-down serial fork of VMEC 8.52 (from STELLOPT `v251`). The goal is to keep only the Stellarator MHD equilibrium solver — nothing else. Preprocessor branches that existed in the upstream code have already been expanded; the source you see is the post-preprocessor code.

The serial-only / educational nature drives several architecture decisions: most state is in module-level globals (no per-thread duplication), and the code is organized to be readable rather than fast.

## Build

CMake out-of-tree build. Required system packages (Debian names): `cmake`, `gfortran`, `libnetcdff-dev`, `libblas-dev`, `liblapack-dev`, `libfftw3-dev`.

```bash
mkdir -p build && cd build
cmake ..
make -j
```

Outputs land in `build/bin/`:
- `xvmec` — main equilibrium solver
- `xnestor` — stand-alone NESTOR (free-boundary vacuum field solver), reads/writes netCDF
- `test_spline_akima` — unit test binary

Submodules `json-fortran` and `abscab-fortran` are required and are wired into the build via `add_subdirectory`. A recursive `git clone` is enough to bring them in; no post-clone cleanup is needed.

## Run

```bash
cd test
../build/bin/xvmec input.solovev      # fixed-boundary smoke test (Solov'ev)
../build/bin/xvmec input.test.vmec    # free-boundary CTH-like (needs mgrid_test.nc)
```

VMEC takes either a full input filename (`input.foo`) or just the extension (`foo`); it extracts `input_extension` from whatever you pass and uses it to name all output files.

Standard outputs in the working directory: `wout_<ext>.nc`, `threed1.<ext>`, `jxbout_<ext>.nc`, `mercier.<ext>`.

`test/coverage/run_all.sh` re-runs the SOLOVEV / HELIOTRON / W7X cases — useful as a broad smoke test after non-trivial changes.

## Tests

```bash
cd build && ctest          # currently runs spline_akima_reflection_symmetry
```

Unit tests live in `test/unit/` and are built as standalone executables that link only the specific sources under test (see `test/unit/CMakeLists.txt`) — they deliberately avoid pulling in the full VMEC library / NetCDF / FFTW / LAPACK stack, so they're a good template for adding new pure-numerics tests.

The bulk of `test/` is reference input files and reference NetCDF outputs (`wout_*.nc`, `jxbout_*.nc`, `mercier.*`) collected from various sources (STELLOPT, DESC, SIMSOPT, …). There is no automated regression harness — comparison is currently manual / via the debug-dump JSON workflow below.

## Architecture: the equilibrium solve

The driver is `src/vmec.f90` (`PROGRAM xvmec` → `subroutine vmec`). The flow is:

1. `reset_params` → `readin` (which calls `read_indata`, `heading`, `read_mgrid`, `allocate_nunv`) → `fixaray` (NS-invariant arrays)
2. **Outer `jacob_off = 0,1` loop** — if the first multi-grid step yields a bad Jacobian, the whole sequence is retried with an inserted `ns=3` initial step.
3. **Multi-grid loop over `ns_array(igrid)`** — each step calls `initialize_radial(nsval, ns_old, …)` (which interpolates the previous solution onto the new mesh) and then **`eqsolve`**, the actual nonlinear equilibrium iteration.
4. `fileout` writes `wout`, `threed1`, `jxbout`, `mercier`.

`eqsolve` (`src/eqsolve.f90`) drives the force-iteration loop, calling **`funct3d`** repeatedly. `funct3d` (`src/funct3d.f90`) is the forward model: given Fourier coefficients of the flux-surface geometry, it evaluates the MHD force residual in Fourier space. This is the inner loop you'll spend most time in when debugging physics.

The free-boundary path inside `funct3d` calls NESTOR (`src/NESTOR/`). NESTOR can run three ways:
- **compiled-in** (default)
- **dump reference I/O** for the internal NESTOR call — set `ldump_vacuum_ref = .true.` (local variable in `funct3d.f90`)
- **out-of-process** — set `lexternal_nestor = .true.` and edit the `nestor_executable` string in `funct3d.f90` to point at `xnestor` or a Python implementation

These three booleans are hard-coded locals in `funct3d.f90`, not input-namelist flags — search for them there when reconfiguring.

## Architecture: state organization

Almost all shared state lives in modules under `src/data/`:
- `vmec_input` — the `INDATA` namelist (everything the user can set, including all the `dump_*` debug flags)
- `vmec_main`, `vmec_params`, `vmec_dim`, `vmec_persistent`, `vparams` — solver state, dimensions, constants
- `xstuff` — the state vectors `xc` (solution) and `scalxc` (scaling)
- `realspace`, `vforces` — real-space fields and forces
- `vmec_io`, `vmercier`, `vsvd0` — output / Mercier / SVD-related state
- `stel_kinds`, `stel_constants` — precision (`rprec`/`dp`) and physical constants

`src/ezcdf/` is the bundled NetCDF wrapper. `src/NESTOR/data/` holds NESTOR-specific module state (`vacmod`, `vacmod0`, `vac_persistent`, `nestor_io`).

## Architecture: the angle constraint

The poloidal angle is not uniquely defined. Upstream STELLOPT had two alternatives switched by the `_HBANGLE` preprocessor flag:
1. Hirshman–Breslau explicit spectrally-optimized series
2. an `m=1`-coefficient constraint (the `lconm1=.true.` path)

In this repo the preprocessor has been pre-expanded with `_HBANGLE` **undefined**, so option 2 (the old `m=1` constraint) is what actually runs. Code regions that originated from the `#ifndef _HBANGLE` branch are still marked with comment fences:

```fortran
! #ifndef _HBANGLE
…
! #end /* ndef _HBANGLE */
```

Treat those fences as documentation, not active preprocessor directives.

## Debugging: the dbgout / JSON workflow

`src/dbgout.f90` provides `open_dbg_context(name, repetition, id)` which writes per-context JSON snapshots when the corresponding `dump_<name>` flag in the `INDATA` namelist is `.true.`. The full list of flags is at `src/data/vmec_input.f90:89` onward. When triggered, files appear in `<input_extension>/<context>/<context>_<ns>_<iteration>_<occurrence>.<ext>.json`.

Each `dump_*` flag corresponds to a string branch in `open_dbg_context` (the big `if/else if` chain in `dbgout.f90`) — adding a new dump point means **both** adding a flag in `vmec_input.f90` and an arm in `dbgout.f90`.

Two global counters guard repetition naming: `vacuum_calls` and `num_eqsolve_retries` (the latter is incremented in exactly one place at the top of `eqsolve.f90` — preserve that invariant).

`test/plot_dbg_json.py` plots quantities/surfaces from these JSON dumps (see README for invocation patterns).

## Conventions worth knowing

- File extensions: `.f90` is free-form, `.f` is fixed-form. Both are present and both are listed explicitly in `src/CMakeLists.txt` — when adding a source file, add it to that list (CMake is not globbing).
- `-cpp` is always on (set in the top-level `CMakeLists.txt`), so `#ifdef NETCDF` etc. still work even though `_HBANGLE` has been expanded out.
- Real precision goes through `rprec` / `dp` from `stel_kinds`; don't hard-code `real(8)`.
- Module `.mod` files are collected under `build/modules/` (and `build/modules/vmec/` for the `vmec` library specifically).
