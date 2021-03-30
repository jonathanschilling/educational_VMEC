cmake_minimum_required(VERSION 3.15.0)

project(vmec LANGUAGES C CXX Fortran)

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

# Fortran specific settings. The first setting tells the compiler to use the C preprocessor. The second places a common directory to place all of the module files.
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

set(SCALAPACK_PREFIX "/usr/local" CACHE PATH "scalapack is installed under this directory")
set(BLACS_PREFIX "/usr/local" CACHE PATH "blacs is installed under this directory")
set(SCALAPACK_LIB_DIR "/usr/local/lib" CACHE PATH "Direct path to where scalapack is installed")
set(BLACS_LIB_DIR "/usr/local/lib" CACHE PATH "Direct path to where blacs is installed")

find_library (SCALAPACK_LIB NAMES scalapack ${SCALAPACK_LIB_NAME} HINTS "${SCALAPACK_LIB_DIR}" "${SCALAPACK_PREFIX}/lib" ENV LD_LIBRARY_PATH)
message(STATUS "scalapack_lib is ${SCALAPACK_LIB}")
if (SCALAPACK_LIB)
# Some builds of SCALAPACK require linking to external BLACS Libraries
    find_library (BLACS_LIB NAMES blacs ${BLACS_LIB_NAME} HINTS "${BLACS_LIB_DIR}" "${BLACS_PREFIX}/lib" ENV LD_LIBRARY_PATH)
	if (BLACS_LIB)
		list (APPEND SCALAPACK_LIB ${BLACS_LIB})
	endif ()
    find_library (BLACSCinit_LIB NAMES blacsCinit ${BLACS_CINIT_NAME} HINTS "${BLACS_LIB_DIR} ${BLACS_PREFIX}/lib" ENV LD_LIBRARY_PATH)
	if (BLACSCinit_LIB)
		list (APPEND SCALAPACK_LIB ${BLACSCinit_LIB})
	endif ()
    find_library (BLACSF77init_LIB NAMES blacsF77init ${BLACS_F77INIT_NAME} HINTS "${BLACS_LIB_DIR} ${BLACS_PREFIX}/lib" ENV LD_LIBRARY_PATH)
	if (BLACSF77init_LIB)
		list (APPEND SCALAPACK_LIB ${BLACSF77init_LIB})
	endif ()
else ()
	#if (USE_PARVMEC) # This option was not used anywhere
    message (FATAL_ERROR "SCALAPACK library is required. A suitable library was not found.")
	#endif ()
endif ()
message(STATUS "updated scalapack_lib is ${SCALAPACK_LIB}")

add_subdirectory(src)
