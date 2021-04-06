
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_attrib.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_GenGet.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_GenPut.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_inqvar.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_opncls.f90
	${CMAKE_CURRENT_SOURCE_DIR}/handle_err.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
