
# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bextrema.f
	${CMAKE_CURRENT_SOURCE_DIR}/blk3d.f
	${CMAKE_CURRENT_SOURCE_DIR}/build_matrices_legendre.f
	${CMAKE_CURRENT_SOURCE_DIR}/CHFDV.f
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
	${CMAKE_CURRENT_SOURCE_DIR}/convert_boundary.f
	${CMAKE_CURRENT_SOURCE_DIR}/convert_boundary_PG.f
	${CMAKE_CURRENT_SOURCE_DIR}/ddsum.f
	${CMAKE_CURRENT_SOURCE_DIR}/erf.f
	${CMAKE_CURRENT_SOURCE_DIR}/factorial.f
	${CMAKE_CURRENT_SOURCE_DIR}/fmax.f
	${CMAKE_CURRENT_SOURCE_DIR}/fmin.f
	${CMAKE_CURRENT_SOURCE_DIR}/getcarg.f
	${CMAKE_CURRENT_SOURCE_DIR}/getfilesize.f
	${CMAKE_CURRENT_SOURCE_DIR}/init_random_seed.f90
	${CMAKE_CURRENT_SOURCE_DIR}/integral.f
	${CMAKE_CURRENT_SOURCE_DIR}/legendre_to_power.f
	${CMAKE_CURRENT_SOURCE_DIR}/line_segment.f
	${CMAKE_CURRENT_SOURCE_DIR}/multxy.f
	${CMAKE_CURRENT_SOURCE_DIR}/parse_extension.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHDF.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHEV.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHEZ.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHFD.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHIM.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHSP.f
	${CMAKE_CURRENT_SOURCE_DIR}/PCHST.f
	${CMAKE_CURRENT_SOURCE_DIR}/point_wise_to_power.f
	${CMAKE_CURRENT_SOURCE_DIR}/power_to_legendre.f
	${CMAKE_CURRENT_SOURCE_DIR}/profile_functions.f
	${CMAKE_CURRENT_SOURCE_DIR}/pxffork.f
	${CMAKE_CURRENT_SOURCE_DIR}/pxfgetpid.f
	${CMAKE_CURRENT_SOURCE_DIR}/pxfwait.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_namelist.f
	${CMAKE_CURRENT_SOURCE_DIR}/second0.f
	${CMAKE_CURRENT_SOURCE_DIR}/solver.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_akima.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_akima_int.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_cubic.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_cubic_int.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline.f
	${CMAKE_CURRENT_SOURCE_DIR}/spline_it.f
	${CMAKE_CURRENT_SOURCE_DIR}/splint.f
	${CMAKE_CURRENT_SOURCE_DIR}/strip_comments.f
	${CMAKE_CURRENT_SOURCE_DIR}/tolower.f90
	${CMAKE_CURRENT_SOURCE_DIR}/tridslv.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_chdir.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_getenv.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_putenv.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_system.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)

