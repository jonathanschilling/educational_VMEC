# Unit tests for individual numerical routines.
#
# Standalone targets that compile only the sources under test, so they do not
# depend on the full vmec library or its NetCDF/FFTW/LAPACK link.

add_executable (test_spline_akima
  ${CMAKE_CURRENT_SOURCE_DIR}/test_spline_akima.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90
  ${CMAKE_SOURCE_DIR}/src/spline_akima.f)
set_target_properties (test_spline_akima PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
add_test (NAME spline_akima_reflection_symmetry COMMAND test_spline_akima)

add_executable (test_cross_section_extrema
  ${CMAKE_CURRENT_SOURCE_DIR}/test_cross_section_extrema.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_constants.f90
  ${CMAKE_SOURCE_DIR}/src/cross_section.f90)
set_target_properties (test_cross_section_extrema PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules_cross_section)
add_test (NAME cross_section_extrema_shift COMMAND test_cross_section_extrema)
