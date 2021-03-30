cmake_minimum_required (VERSION 3.15.0)

#message(STATUS "CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

# Start Project and note the language used. The BLAS and LAPACK libraries for ifort require a working C and C++ compiler.
project (VMEC Fortran C CXX)

#enable_testing ()


# Activate or Deactivate code profiling. The profiler maybe activated by setting -DUSE_PROFILER=ON on the cmake command line.
set (vmec_sources "")
# List all the subdirectories these are in library dependancy order.
#add_subdirectory (Sources)

# Add subdirectories for all the sources.
add_subdirectory (LIBSTELL)
add_subdirectory (General)
add_subdirectory (Hessian)
add_subdirectory (Initialization_Cleanup)
add_subdirectory (Input_Output)
add_subdirectory (NESTOR_vacuum)
add_subdirectory (Splines)
add_subdirectory (TimeStep)

set (vmec_sources "${vmec_sources}")# PARENT_SCOPE)

add_library (vmec STATIC ${vmec_sources})
set_target_properties (vmec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
target_include_directories(vmec PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/TimeStep/vmec.f)
target_link_libraries (xvmec vmec)

target_link_libraries (vmec PUBLIC ${NETCDF_LIBRARIES})

if (SCALAPACK_LIB)
	target_link_libraries (vmec PUBLIC ${SCALAPACK_LIB})
endif ()
if (BLAS_FOUND)
	target_link_libraries (vmec PUBLIC ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (vmec PUBLIC ${LAPACK_LIBRARIES})
endif ()

get_target_property(VMEC_LINK_LIB vmec LINK_LIBRARIES)
message(STATUS "vmec linked libraries are ${VMEC_LINK_LIB}")
set(VMEC_LINK_LIB ${VMEC_LINK_LIB} PARENT_SCOPE)

