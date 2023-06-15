## Educational VMEC

This is a heavily stripped-down version of the serial implementation of VMEC 8.52.
It is forked from the `v251` branch of the [STELLOPT](https://github.com/PrincetonUniversity/STELLOPT) repository.

The goal of this project is to have a version of VMEC
which only computes the Stellarator MHD equilibrium and nothing more.

The `cmake` build system for stand-alone VMEC is borrowed
from [hiddenSymmetries/VMEC2000](https://github.com/hiddenSymmetries/VMEC2000)
and from [ORNL-Fusion/LIBSTELL](https://github.com/ORNL-Fusion/LIBSTELL).

## Building
Required packages on Debian stable are:
`cmake`, `gfortran`, `libnetcdff-dev`, `libblas-dev`, `liblapack-dev`, `libfftw3-dev`

You need to download this repository and the included `json-fortran`.
The easiest way to do this is via a recursive `git clone`:

```bash
> git clone --recursive https://github.com/jonathanschilling/educational_VMEC.git
```

Ignore warnings about missing access to `src/vac2` and `src/vac3`.
After cloning, the `json-fortran` directory is somewhat messed up and you need to fix that:

```bash
> cd educational_VMEC/json-fortran
> git restore --staged .
> git checkout .
> cd ..
```

This is a fairly standard CMake setup, if you are used to it.
Here is how it works:
 * Create a directory `build` in the main folder: `mkdir build`
 * Go into the `build` directory: `cd build`
 * Run CMake: `cmake ..`
 * Execute the actual build process: `make` (optional multi-threaded build: `make -j`)
 * The VMEC executable `xvmec` is then located in `build/bin` with respect to the main folder.
 
## Example Execution
 * Change into the `test` dir: `cd test`
 * Run the [Solov'ev test case](https://princetonuniversity.github.io/FOCUS/notes/Coil_design_codes_benchmark.html#Equiblirium--): `../build/bin/xvmec input.solovev`

## Debug Output

The [`dbgout`](src/dbgout.f90) module allows to (optionally) write a bunch of data from within the code into separate JSON files.
Writing these additional quantities is disabled by default and can be enabled by logical flags in the `INDATA` namelist.
A full list of the available flags (ordered by module where they are used) is available in [vmec_input.f90](src/data/vmec_input.f90#L89).
A folder named the same as the extension of the input file will be created.
In there, subfolders named after each of the enabled `dump_*` flags will be created (`dump_forces=.true.` --> `forces/`)
In each of these subfolders, JSON files will be created according to the following naming scheme:
`forces_<ns>_<iteration>_<occurence>.<extension>.json` where
- `ns` is the number of flux surfaces of the respective multi-grid step,
- `iteration` is the iteration number printed to the screen and
- `occurence` is a linear counter used in case the respective quantities need to get dumped more than once in the respective context (e.g. before and after some important subroutine call)

As a first start, I usually look at JSON files using the Firefox browser, which has a rather nice built-in JSON viewer.
For plotting, I have written a Python utility [`plot_dbg_json.py`](test/plot_dbg_json.py),
which can plot the quantities from the JSON files produced by `educational_VMEC`.

Use it as follows:
1. to plot all quantities (at all flux surfaces) in a given JSON file:
 
   ```
   ../../plot_dbg_json.py forces_00015_000001_01.test.vmec.json
   ```
   
2. the list of quantities to plot can be specified as follows (for quicker iterations, since plotting everything can take a while):

    ```
    ../../plot_dbg_json.py forces_00015_000001_01.test.vmec.json --quantities ru12
    ../../plot_dbg_json.py forces_00015_000001_01.test.vmec.json --quantities ru12 zu12
    ```
   
3. to plot a subset of flux surfaces (0-based indexing):
 
   ```
   ../../plot_dbg_json.py forces_00015_000001_01.test.vmec.json --surfaces 6
   ../../plot_dbg_json.py forces_00015_000001_01.test.vmec.json --surfaces 0 5 14
   ```
   
4. to plot a subset of quantities at a subset of surfaces:

   ```
   ../../plot_dbg_json.py forces_00015_000001_01.test.vmec.json --quantities ru12 zu12 --surfaces 6 7
   ```

## External NESTOR
The free-boundary part of VMEC is the Neumann Solver for Toroidal Systems (NESTOR).
Its source code is in a separate folder [`NESTOR`](src/NESTOR).
The appropriate reference is https://doi.org/10.1016/0021-9991(86)90055-0 .

This version of NESTOR can be run stand-alone. It reads its inputs from a netCDF file and writes its outputs into another netCDF file.
The main executable of this stand-alone version of NESTOR is [`nestor_main.f90`](src/NESTOR/nestor_main.f90).
The input and output files are read and written in [`nestor_io.f90`](src/NESTOR/data/nestor_io.f90).

This version of VMEC can be configured to dump the corresponding input and output files, but still run the compiled-in version of NESTOR.
This is enabled via the logical flag `ldump_vacuum_ref` in [`funct3d.f90`](src/funct3d.f90).

Also, an external NESTOR implementation can be called instead of using the compiled-in version of NESTOR.
This is enabled via the logical flag `lexternal_nestor` in [`funct3d.f90`](src/funct3d.f90).
The corresponding system call to execute the external NESTOR implementation has to be specified in
`nestor_executable` in [`funct3d.f90`](src/funct3d.f90).

## Angle Constraint
The poloidal angle-like coordinate is a priori not uniquely defined and needs special care.
The version of VMEC from the STELLOPT repo had essentially two options for this.
They were alternatively compiled in via the preprocessor flag `_HBANGLE`.

1. The Hirshman-Breslau explicit spectrally optimized Fourier series (see https://doi.org/10.1063/1.872954 for details) and
2. an unknown mixture of several constraints of the `m=1` Fourier coefficients (the logical `lconm1` is true for this constraint).

By default, the `_HBANGLE` preprocessor flag is not active and thus, the "old" `m=1` constraint is active.

This version of VMEC has all preprocessor flags already expanded.
It became clear that it is nevertheless useful to have at least a vague idea of what parts of the code are related to the angle constraint.
Therefore, those parts of VMEC related to the `m=1`constraint are marked to start with 

```Fortran
! #ifndef _HBANGLE
```

and end with


```Fortran
! #end /* ndef _HBANGLE */
```

Note that the [original documentation](https://github.com/jonathanschilling/educational_VMEC/blob/master/vmec_info.md) states:
> THE POLOIDAL ANGLE IS DETERMINED BY MINIMIZING <M> = m\*\*2 S(m) , WHERE S(m) = Rm\*\*2 + Zm\*\*2 .
