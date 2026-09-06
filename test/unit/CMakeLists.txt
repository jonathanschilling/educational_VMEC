# Numerical unit tests and a solver-level initialization regression.
#
# The spline target compiles only the sources under test, so it does not
# depend on the full vmec library or its NetCDF/FFTW/LAPACK link.

add_executable (test_spline_akima
  ${CMAKE_CURRENT_SOURCE_DIR}/test_spline_akima.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90
  ${CMAKE_SOURCE_DIR}/src/spline_akima.f)
set_target_properties (test_spline_akima PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
add_test (NAME spline_akima_reflection_symmetry COMMAND test_spline_akima)

add_executable(check_jacobian_recovery
  check_jacobian_recovery.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90)
set_target_properties(check_jacobian_recovery PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/recovery_modules)
target_link_libraries(check_jacobian_recovery PRIVATE NetCDF::NetCDF)
add_test(NAME initial_jacobian_recovery
  COMMAND ${CMAKE_COMMAND}
    -DVMEC=$<TARGET_FILE:xvmec>
    -DCHECKER=$<TARGET_FILE:check_jacobian_recovery>
    -DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/input.initial_jacobian_recovery
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/jacobian_recovery
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_jacobian_recovery.cmake)
