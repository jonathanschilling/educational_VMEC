cmake_minimum_required(VERSION 3.15.0)

project(vmec LANGUAGES C CXX Fortran)

# Start Project and note the language used. The BLAS and LAPACK libraries for ifort require a working C and C++ compiler.
project (VMEC Fortran C CXX)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package (NetCDF MODULE REQUIRED)

find_package (FFTW MODULE REQUIRED COMPONENTS DOUBLE_LIB)

set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

# Set a directories to build all binary files.
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)
set (LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Add extra compiler flags to various builds.
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_RELEASE_FLAGS}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_DEBUG_FLAGS}")
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FLAGS}")

# Fortran specific settings. The first setting tells the compiler to use the C preprocessor.
# The second places a common directory to place all of the module files.
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

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

if (FFTW_DOUBLE_LIB_FOUND)
	message (STATUS "Using fftw")
	include_directories (${FFTW_INCLUDE_DIRS})
else ()
	message (FATAL_ERROR "Failed to find the required fftw library.")
endif ()

set (vmec_sources "")
add_subdirectory(src)
add_subdirectory(json-fortran)
add_subdirectory(abscab-fortran)

# if you have access to the vac2 and vac3 versions of NESTOR (contact me), they are included here
#add_subdirectory(src/vac2)
#add_subdirectory(src/vac3)

add_library (vmec STATIC ${vmec_sources})
set_target_properties (vmec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)
target_include_directories(vmec PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}/vmec)

target_link_libraries (vmec PUBLIC ${NETCDF_LIBRARIES})
target_link_libraries (vmec PUBLIC ${FFTW_DOUBLE_LIB})

if (BLAS_FOUND)
	target_link_libraries (vmec PUBLIC ${BLAS_LIBRARIES})
endif ()
if (LAPACK_FOUND)
	target_link_libraries (vmec PUBLIC ${LAPACK_LIBRARIES})
endif ()

get_target_property(VMEC_LINK_LIB vmec LINK_LIBRARIES)
message(STATUS "vmec linked libraries are ${VMEC_LINK_LIB}")



# Define an executable and link all libraries.
add_executable (xvmec ${CMAKE_CURRENT_SOURCE_DIR}/src/vmec.f90)
target_link_libraries (xvmec vmec)

### stand-alone NESTOR
add_executable (xnestor ${CMAKE_CURRENT_SOURCE_DIR}/src/NESTOR/nestor_main.f90)
target_link_libraries (xnestor vmec)

