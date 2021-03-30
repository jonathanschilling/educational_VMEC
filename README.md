## Educational VMEC

This is a heavily stripped-down version of the serial implementation of VMEC 8.52.
It is forked from the `v251` branch of the [STELLOPT](https://github.com/PrincetonUniversity/STELLOPT) repository.

The goal of this project is to have a version of VMEC
which only computes the Stellarator MHD equilibrium and nothing more.

The `cmake` build system for stand-alone VMEC is borrowed
from [hiddenSymmetries/VMEC2000](https://github.com/hiddenSymmetries/VMEC2000)
and from [ORNL-Fusion/LIBSTELL](https://github.com/ORNL-Fusion/LIBSTELL).

## Building
This is a fairly standard CMake setup, if you are used to it.
Here is how it works:
 * Create a directory `build` in the main folder: `mkdir build`
 * Go into the `build` directory: `cd build`
 * Run CMake: `cmake ..`
 * Execute the actual build process: `make` (optional multi-threaded build: `make -j`)
 * The VMEC executable `xvmec` is then located in `build/build/bin` with respect to the main folder.

## Example Execution
 * Change into the `test` dir: `cd test`
 * Run the [Solov'ev test case](https://princetonuniversity.github.io/FOCUS/notes/Coil_design_codes_benchmark.html#Equiblirium--): `../build/build/bin/xvmec input.solovev`
