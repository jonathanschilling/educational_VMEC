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
 * The VMEC executable `xvmec` is then located in `build/bin` with respect to the main folder.
 
## Example Execution
 * Change into the `test` dir: `cd test`
 * Run the [Solov'ev test case](https://princetonuniversity.github.io/FOCUS/notes/Coil_design_codes_benchmark.html#Equiblirium--): `../build/bin/xvmec input.solovev`

## Angle Constraint
The poloidal angle-like coordinate is a priori not uniquely defined and needs special care.
The original VMEC from the STELLOPT repo had essentially two options for this.
They were alternatively compiled in via the preprocessor flag `_HBANGLE`.
1. The Hirshman-Breslau explicit spectrally optimized Fourier series (see https://doi.org/10.1063/1.872954 for details) and
2. an unknown mixture of several constraints of the `m=1` Fourier coefficients (the logical `lconm1` is true for this constraint).
By default, the `_HBANGLE` preprocessor flag is not active and thus, the "old" `m=1` constraint is active.

This version of VMEC has most, if not all, of its preprocessor flags explicitly expanded.
It became clear that it is nevertheless useful to have at least a vague idea of what parts of the code are related to the angle constraint.
Therefore, those parts of VMEC related to the `m=1`constraint are marked to start with 

```Fortran
! # else ifndef _HBANGLE
```

and end with


```Fortran
! # end ifndef _HBANGLE
```
