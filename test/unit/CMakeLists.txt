# Unit tests for individual numerical routines.
#
# test_spline_akima compiles only the sources under test, so it does not
# depend on the full vmec library or its NetCDF/FFTW/LAPACK link.

add_executable (test_spline_akima
  ${CMAKE_CURRENT_SOURCE_DIR}/test_spline_akima.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90
  ${CMAKE_SOURCE_DIR}/src/spline_akima.f)
set_target_properties (test_spline_akima PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
add_test (NAME spline_akima_reflection_symmetry COMMAND test_spline_akima)

# NESTOR on two current rings inside a circular torus, whose total vacuum field
# on the boundary is zero, on the axisymmetric grid, on three-dimensional grids
# of even and odd nzeta, and on the full poloidal grid of lasym.
add_executable (test_nestor_two_rings
  ${CMAKE_CURRENT_SOURCE_DIR}/test_nestor_two_rings.f90)
target_link_libraries (test_nestor_two_rings PRIVATE vmec)
add_test (NAME nestor_two_rings_axisymmetric COMMAND test_nestor_two_rings 1 0 1 22 F)
add_test (NAME nestor_two_rings_even_nzeta COMMAND test_nestor_two_rings 2 4 32 22 F)
add_test (NAME nestor_two_rings_odd_nzeta COMMAND test_nestor_two_rings 2 4 31 22 F)
add_test (NAME nestor_two_rings_lasym COMMAND test_nestor_two_rings 2 4 31 22 T)

# The same problem on finer grids, where the error that remains is the
# third-order defect of the sums of greenf at the source point.
add_test (NAME nestor_two_rings_source_point_axisymmetric COMMAND test_nestor_two_rings 1 0 1 32 F 6e-5)
add_test (NAME nestor_two_rings_source_point_even_nzeta COMMAND test_nestor_two_rings 2 4 32 44 F 3e-5)
add_test (NAME nestor_two_rings_source_point_odd_nzeta COMMAND test_nestor_two_rings 2 4 31 44 F 3e-5)
add_test (NAME nestor_two_rings_source_point_lasym COMMAND test_nestor_two_rings 2 4 32 44 T 3e-5)

# The axisymmetric grid with nfp other than 1, whose integral in v goes around
# the whole torus.
add_test (NAME nestor_two_rings_axisymmetric_nfp2 COMMAND test_nestor_two_rings 2 0 1 22 F)
add_test (NAME nestor_two_rings_axisymmetric_nfp5 COMMAND test_nestor_two_rings 5 0 1 22 F)
add_test (NAME nestor_two_rings_axisymmetric_nfp5_lasym COMMAND test_nestor_two_rings 5 0 1 22 T)
