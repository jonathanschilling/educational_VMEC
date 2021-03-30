
# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/drive_gmres.f
	${CMAKE_CURRENT_SOURCE_DIR}/gmres.f
	${CMAKE_CURRENT_SOURCE_DIR}/gmresr.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_gmres.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
