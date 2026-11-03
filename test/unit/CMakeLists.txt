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

# A lasym run and the same run shifted vertically must give the same equilibrium.
add_executable (check_vertical_shift
  ${CMAKE_CURRENT_SOURCE_DIR}/check_vertical_shift.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90)
set_target_properties (check_vertical_shift PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/check_vertical_shift_modules)
target_link_libraries (check_vertical_shift PRIVATE NetCDF::NetCDF)
add_test (NAME lasym_vertical_shift
  COMMAND ${CMAKE_COMMAND}
    -DVMEC=$<TARGET_FILE:xvmec>
    -DCHECKER=$<TARGET_FILE:check_vertical_shift>
    -DINPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/vertical_shift
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_vertical_shift.cmake)
