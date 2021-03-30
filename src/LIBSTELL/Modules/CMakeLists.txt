
# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/bcoils_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/biotsavart.f
	${CMAKE_CURRENT_SOURCE_DIR}/boundary_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/bsc_cdf.f
	${CMAKE_CURRENT_SOURCE_DIR}/bsc_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
	${CMAKE_CURRENT_SOURCE_DIR}/coordinate_utilities.f
	${CMAKE_CURRENT_SOURCE_DIR}/cyl_flux.f
	${CMAKE_CURRENT_SOURCE_DIR}/date_and_computer.f
	${CMAKE_CURRENT_SOURCE_DIR}/extcurz_T.f
	${CMAKE_CURRENT_SOURCE_DIR}/fftpack.f
	${CMAKE_CURRENT_SOURCE_DIR}/flin1_z.f
	${CMAKE_CURRENT_SOURCE_DIR}/flint_z.f
	${CMAKE_CURRENT_SOURCE_DIR}/functions.f
	${CMAKE_CURRENT_SOURCE_DIR}/linear_interp.f90
	${CMAKE_CURRENT_SOURCE_DIR}/math_utilities.f
	${CMAKE_CURRENT_SOURCE_DIR}/mgrid_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/modular_coils.f
	${CMAKE_CURRENT_SOURCE_DIR}/profiler.f
	${CMAKE_CURRENT_SOURCE_DIR}/read_wout_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/safe_open_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/stel_constants.f
	${CMAKE_CURRENT_SOURCE_DIR}/stel_kinds.f
	${CMAKE_CURRENT_SOURCE_DIR}/system_mod.f
	${CMAKE_CURRENT_SOURCE_DIR}/v3f_vmec_comm.f
	${CMAKE_CURRENT_SOURCE_DIR}/vcoilpts.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_input.f
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_utils.f
	${CMAKE_CURRENT_SOURCE_DIR}/vparams.f
	${CMAKE_CURRENT_SOURCE_DIR}/vsvd0.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)

