# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/eqsolve.f
	${CMAKE_CURRENT_SOURCE_DIR}/evolve.f
	${CMAKE_CURRENT_SOURCE_DIR}/interp.f
	${CMAKE_CURRENT_SOURCE_DIR}/restart_iter.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
