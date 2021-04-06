	# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bss.f90
	${CMAKE_CURRENT_SOURCE_DIR}/elongation.f90
	${CMAKE_CURRENT_SOURCE_DIR}/eqfor.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fileout.f90
	${CMAKE_CURRENT_SOURCE_DIR}/freeb_data.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fsym_fft.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fsym_invfft.f90
	${CMAKE_CURRENT_SOURCE_DIR}/getbsubs.f90
	${CMAKE_CURRENT_SOURCE_DIR}/getcurmid.f90
	${CMAKE_CURRENT_SOURCE_DIR}/jxbforce.f90
	${CMAKE_CURRENT_SOURCE_DIR}/mercier.f90
	${CMAKE_CURRENT_SOURCE_DIR}/printout.f90
	${CMAKE_CURRENT_SOURCE_DIR}/readin.f90
	${CMAKE_CURRENT_SOURCE_DIR}/read_indata.f90
	${CMAKE_CURRENT_SOURCE_DIR}/timer_sub.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_io.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmercier.f90
	${CMAKE_CURRENT_SOURCE_DIR}/wrout.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)
