
# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/cfode.f
	${CMAKE_CURRENT_SOURCE_DIR}/dlsode_aux1.f
	${CMAKE_CURRENT_SOURCE_DIR}/dlsode_aux2.f
	${CMAKE_CURRENT_SOURCE_DIR}/dlsode.f
	${CMAKE_CURRENT_SOURCE_DIR}/ewset.f
	${CMAKE_CURRENT_SOURCE_DIR}/intdy.f
	${CMAKE_CURRENT_SOURCE_DIR}/lsode.f
	${CMAKE_CURRENT_SOURCE_DIR}/prepj.f
	${CMAKE_CURRENT_SOURCE_DIR}/solsy.f
	${CMAKE_CURRENT_SOURCE_DIR}/stode.f
	${CMAKE_CURRENT_SOURCE_DIR}/vnorm.f
	${CMAKE_CURRENT_SOURCE_DIR}/xerrwv.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
