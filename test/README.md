# Collection of Input Files for VMEC

This is a collection of input files for VMEC that are freely available from various repositories.

## Fixed-Boundary

* NCSX (from https://princetonuniversity.github.io/STELLOPT/VMEC%20Fixed%20Boundary%20Run)

## Free-Boundary

* NCSX (from https://princetonuniversity.github.io/STELLOPT/VMEC%20Free%20Boundary%20Run)

## CTH-like
VMEC run of input.test.vmec (CTH-like free-boundary stellarator run)

1. ldump_vacuum_ref = .true., lexternal_nestor = .false.
--> vac_ref_test.vmec.tar (folder vac_ref as in NESTOR README.md)

2. ldump_vacuum_ref = .false., lexternal_nestor = .true.,
   nestor_executable = "/data2/jonathan/work/code/educational_VMEC/build/bin/xnestor"
   (original Fortran implementation, but run stand-alone)
--> vac_internal_test.vmec.tar (folder vac as in NESTOR README.md)

3. ldump_vacuum_ref = .false., lexternal_nestor = .true.,
   nestor_executable = "python3 /home/IPP-HGW/jons/work/code/NESTOR/src/main/python/NESTOR.py"
   (Python implementation of NESTOR)
--> vac_test.vmec.tar (folder vac as in NESTOR README.md)
