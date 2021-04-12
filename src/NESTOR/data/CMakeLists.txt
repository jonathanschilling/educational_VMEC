
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/vacmod.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vacmod0.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vac_persistent.f90
	${CMAKE_CURRENT_SOURCE_DIR}/nestor_io.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
