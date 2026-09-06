`input.initial_jacobian_recovery` imposes the outermost boundary of
[QUASR configuration 14802](https://quasr.flatironinstitute.org/simsopt_serials/0014/serial0014802.json)
at `mpol=ntor=6`, with zero pressure and prescribed zero toroidal current.
The toroidal flux comes from the prescribed coil vector potential. Its boundary
coefficients were converted with SIMSOPT and are not changed during recovery.

The cold `[8,16,31]` sequence has a bad initial Jacobian after the axis search.
An explicit `[3,8,16,31]` sequence converges, using `ftol=1e-4` only for the
three-surface bootstrap and `1e-9` for the requested grids. The regression test
requires the automatic recovery to reach the requested `ns=31` and tolerance.
