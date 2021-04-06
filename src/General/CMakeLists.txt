# Append all source files to variable stell_sources. As new files are added this must be updated.
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/add_fluxes.f90
	${CMAKE_CURRENT_SOURCE_DIR}/alias.f90
	${CMAKE_CURRENT_SOURCE_DIR}/aspectratio.f90
	${CMAKE_CURRENT_SOURCE_DIR}/bcovar.f90
	${CMAKE_CURRENT_SOURCE_DIR}/convert.f90
	${CMAKE_CURRENT_SOURCE_DIR}/fbal.f90
	${CMAKE_CURRENT_SOURCE_DIR}/forces.f90
	${CMAKE_CURRENT_SOURCE_DIR}/funct3d.f90
	${CMAKE_CURRENT_SOURCE_DIR}/getfsq.f90
	${CMAKE_CURRENT_SOURCE_DIR}/jacobian.f90
	${CMAKE_CURRENT_SOURCE_DIR}/lamcal.f90
	${CMAKE_CURRENT_SOURCE_DIR}/precondn.f90
	${CMAKE_CURRENT_SOURCE_DIR}/realspace.f90
	${CMAKE_CURRENT_SOURCE_DIR}/residue.f90
	${CMAKE_CURRENT_SOURCE_DIR}/scalfor.f90
	${CMAKE_CURRENT_SOURCE_DIR}/spectrum.f90
	${CMAKE_CURRENT_SOURCE_DIR}/symforce.f90
	${CMAKE_CURRENT_SOURCE_DIR}/symrzl.f90
	${CMAKE_CURRENT_SOURCE_DIR}/tomnsp_mod.f90
	${CMAKE_CURRENT_SOURCE_DIR}/totzsp_mod.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vforces.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_dim.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_main.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_params.f90
	${CMAKE_CURRENT_SOURCE_DIR}/vmec_persistent.f90
	${CMAKE_CURRENT_SOURCE_DIR}/xstuff.f90
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)

set_property(SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/fbal.f APPEND PROPERTY OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/vforce.f")
