cmake_minimum_required(VERSION 3.15.0)

project(vmec LANGUAGES C CXX Fortran)

# Start Project and note the language used. The BLAS and LAPACK libraries for ifort require a working C and C++ compiler.
project (VMEC Fortran C CXX)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package (NetCDF MODULE REQUIRED)

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/build/lib)

# Add extra compiler flags to various builds.
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_RELEASE_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_DEBUG_FLAGS}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS}")

# Fortran specific settings. The first setting tells the compiler to use the C preprocessor.
# The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/build/modules)

find_package (BLAS REQUIRED)
if (NOT ${BLA_VENDOR} MATCHES "^Intel")
    find_package (LAPACK REQUIRED)
endif()

if (NetCDF_FOUND)
	message (STATUS "Using netcdf")
	add_definitions (-DNETCDF)
	include_directories (${NetCDF_INCLUDE_DIRS})

	list (APPEND NETCDF_LIBRARIES "${NetCDF_LIBRARIES}")
else ()
	message (FATAL_ERROR "Failed to find the required netcdf libraries.")
endif ()

set (vmec_sources "")
add_subdirectory(src)

add_library (vmec STATIC ${vmec_sources})
set_target_properties (vmec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
target_include_directories(vmec PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)

# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/src/vmec.f90)
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

