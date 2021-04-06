# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_funct3d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_ns.f90
	${CMAKE_CURRENT_SOURCE_DIR}/allocate_nunv.f90
	${CMAKE_CURRENT_SOURCE_DIR}/close_all_files.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fixaray.f90
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_funct3d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_ns.f90
	${CMAKE_CURRENT_SOURCE_DIR}/free_mem_nunv.f90
	${CMAKE_CURRENT_SOURCE_DIR}/free_persistent_mem.f90
	${CMAKE_CURRENT_SOURCE_DIR}/guess_axis.f90
	${CMAKE_CURRENT_SOURCE_DIR}/heading.f90
	${CMAKE_CURRENT_SOURCE_DIR}/init_geometry.f90
	${CMAKE_CURRENT_SOURCE_DIR}/initialize_radial.f90
	${CMAKE_CURRENT_SOURCE_DIR}/magnetic_fluxes.f90
	${CMAKE_CURRENT_SOURCE_DIR}/open_output_files.f90
	${CMAKE_CURRENT_SOURCE_DIR}/profil1d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/profil3d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/reset_params.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vsetup.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
