
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/realspace.f90
	${CMAKE_CURRENT_SOURCE_DIR}/stel_constants.f90
	${CMAKE_CURRENT_SOURCE_DIR}/stel_kinds.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vforces.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_dim.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_input.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_io.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_main.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_params.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_persistent.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmercier.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vparams.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vsvd0.f90
	${CMAKE_CURRENT_SOURCE_DIR}/xstuff.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
