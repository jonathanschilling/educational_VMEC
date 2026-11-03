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

# A vertically stable lasym ring started above its midplane must descend.
add_executable (check_raised_ring
  ${CMAKE_CURRENT_SOURCE_DIR}/check_raised_ring.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90)
set_target_properties (check_raised_ring PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/check_raised_ring_modules)
target_link_libraries (check_raised_ring PRIVATE NetCDF::NetCDF)
add_executable (make_raised_ring_mgrid
  ${CMAKE_CURRENT_SOURCE_DIR}/make_raised_ring_mgrid.f90
  ${CMAKE_SOURCE_DIR}/src/data/stel_kinds.f90)
set_target_properties (make_raised_ring_mgrid PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/make_raised_ring_mgrid_modules)
target_link_libraries (make_raised_ring_mgrid PRIVATE NetCDF::NetCDF)
add_test (NAME lasym_raised_ring
  COMMAND ${CMAKE_COMMAND}
    -DVMEC=$<TARGET_FILE:xvmec>
    -DCHECKER=$<TARGET_FILE:check_raised_ring>
    -DMGRID_WRITER=$<TARGET_FILE:make_raised_ring_mgrid>
    -DINPUT_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/raised_ring
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_raised_ring.cmake)
