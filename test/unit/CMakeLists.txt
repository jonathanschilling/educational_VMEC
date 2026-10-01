# Unit tests for individual numerical routines, and a solver-level regression.
#
# The unit tests are standalone targets that compile only the sources under
# test, so they do not depend on the full vmec library or its
# NetCDF/FFTW/LAPACK link.

add_executable (test_spline_akima
  ${CMAKE_CURRENT_SOURCE_DIR}/test_spline_akima.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90
  ${CMAKE_SOURCE_DIR}/src/spline_akima.f)
set_target_properties (test_spline_akima PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
add_test (NAME spline_akima_reflection_symmetry COMMAND test_spline_akima)

# A fixed-boundary run that restarts after bad Jacobians and must converge.
add_executable (check_converged
  ${CMAKE_CURRENT_SOURCE_DIR}/check_converged.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90)
set_target_properties (check_converged PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/check_converged_modules)
target_link_libraries (check_converged PRIVATE NetCDF::NetCDF)
add_test (NAME evaluated_backup_restart
  COMMAND ${CMAKE_COMMAND}
    -DVMEC=$<TARGET_FILE:xvmec>
    -DCHECKER=$<TARGET_FILE:check_converged>
    -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/input.evaluated_backup
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/evaluated_backup
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_converges.cmake)
