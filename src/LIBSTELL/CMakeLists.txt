
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_attrib.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_GenGet.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_GenPut.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_inqvar.f90
	${CMAKE_CURRENT_SOURCE_DIR}/ezcdf_opncls.f90
	${CMAKE_CURRENT_SOURCE_DIR}/handle_err.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/bextrema.f
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
	${CMAKE_CURRENT_SOURCE_DIR}/fmax.f
	${CMAKE_CURRENT_SOURCE_DIR}/fmin.f
	${CMAKE_CURRENT_SOURCE_DIR}/getcarg.f
	${CMAKE_CURRENT_SOURCE_DIR}/line_segment.f
	${CMAKE_CURRENT_SOURCE_DIR}/parse_extension.f
	${CMAKE_CURRENT_SOURCE_DIR}/profile_functions.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_namelist.f
	${CMAKE_CURRENT_SOURCE_DIR}/second0.f
	${CMAKE_CURRENT_SOURCE_DIR}/solver.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_akima.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_akima_int.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_cubic.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_cubic_int.f
	${CMAKE_CURRENT_SOURCE_DIR}/tolower.f90
	${CMAKE_CURRENT_SOURCE_DIR}/tridslv.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_getenv.f
	${CMAKE_CURRENT_SOURCE_DIR}/biotsavart.f
	${CMAKE_CURRENT_SOURCE_DIR}/bsc_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
	${CMAKE_CURRENT_SOURCE_DIR}/cyl_flux.f
	${CMAKE_CURRENT_SOURCE_DIR}/date_and_computer.f
 	${CMAKE_CURRENT_SOURCE_DIR}/functions.f
	${CMAKE_CURRENT_SOURCE_DIR}/mgrid_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_wout_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/safe_open_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/stel_constants.f
	${CMAKE_CURRENT_SOURCE_DIR}/stel_kinds.f
	${CMAKE_CURRENT_SOURCE_DIR}/system_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_utils.f
	${CMAKE_CURRENT_SOURCE_DIR}/vparams.f
	${CMAKE_CURRENT_SOURCE_DIR}/vsvd0.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)

