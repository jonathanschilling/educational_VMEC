module dbgout
  use json
  implicit none

contains

function open_dbg_context(context_name)
  use vmec_dim,   only: ns
  use vmec_main,  only: iter2
  use vmec_input
  implicit none
  
  character(len=*), intent(in) :: context_name
  logical :: open_dbg_context
  
  character(len=255) :: dump_filename
  logical            :: should_write
  
  ! check if requested context is enabled by input flags
  if      (trim(context_name) .eq. "add_fluxes") then
    open_dbg_context         = dump_add_fluxes
  else if (trim(context_name) .eq. "metric") then
    open_dbg_context         = dump_metric
  else if (trim(context_name) .eq. "volume") then
    open_dbg_context         = dump_volume
  else if (trim(context_name) .eq. "bcontrav") then
    open_dbg_context         = dump_bcontrav
  else if (trim(context_name) .eq. "bcov") then
    open_dbg_context         = dump_bcov
  else if (trim(context_name) .eq. "lambda_forces") then
    open_dbg_context         = dump_lambda_forces
  else if (trim(context_name) .eq. "bcov_full") then
    open_dbg_context         = dump_bcov_full
  else if (trim(context_name) .eq. "precondn") then
    open_dbg_context         = dump_precondn
  else if (trim(context_name) .eq. "forceNorms_tcon") then
    open_dbg_context         = dump_forceNorms_tcon
  else if (trim(context_name) .eq. "lulv_comb") then
    open_dbg_context         = dump_lulv_comb
  else if (trim(context_name) .eq. "calc_fbal") then
    open_dbg_context         = dump_calc_fbal
  else if (trim(context_name) .eq. "evolve") then
    open_dbg_context         = dump_evolve
  else if (trim(context_name) .eq. "fixaray") then
    open_dbg_context         = dump_fixaray
  else if (trim(context_name) .eq. "spectral_constraint") then
    open_dbg_context         = dump_spectral_constraint
  else if (trim(context_name) .eq. "forces") then
    open_dbg_context         = dump_forces
  else if (trim(context_name) .eq. "funct3d_geometry") then
    open_dbg_context         = dump_funct3d_geometry
  else if (trim(context_name) .eq. "rbsq") then
    open_dbg_context         = dump_rbsq
  else if (trim(context_name) .eq. "constraint_force") then
    open_dbg_context         = dump_constraint_force
  else if (trim(context_name) .eq. "guess_axis") then
    open_dbg_context         = dump_guess_axis
  else if (trim(context_name) .eq. "interp") then
    open_dbg_context         = dump_interp
  else if (trim(context_name) .eq. "jacobian") then
    open_dbg_context         = dump_jacobian
  else if (trim(context_name) .eq. "lamcal") then
    open_dbg_context         = dump_lamcal
  else if (trim(context_name) .eq. "profil1d") then
    open_dbg_context         = dump_profil1d
  else if (trim(context_name) .eq. "profil3d") then
    open_dbg_context         = dump_profil3d
  else if (trim(context_name) .eq. "readin_boundary") then
    open_dbg_context         = dump_readin_boundary
  else if (trim(context_name) .eq. "fsq") then
    open_dbg_context         = dump_fsq
  else if (trim(context_name) .eq. "scale_m1") then
    open_dbg_context         = dump_scale_m1
  else if (trim(context_name) .eq. "scalfor_out") then
    open_dbg_context         = dump_scalfor_out
  else if (trim(context_name) .eq. "fsq1") then
    open_dbg_context         = dump_fsq1
  else if (trim(context_name) .eq. "scalfor_R") then
    open_dbg_context         = dump_scalfor_R
  else if (trim(context_name) .eq. "scalfor_Z") then
    open_dbg_context         = dump_scalfor_Z
  else if (trim(context_name) .eq. "symforce") then
    open_dbg_context         = dump_symforce
  else if (trim(context_name) .eq. "tomnsps") then
    open_dbg_context         = dump_tomnsps
  else if (trim(context_name) .eq. "tomnspa") then
    open_dbg_context         = dump_tomnspa
  else if (trim(context_name) .eq. "multigrid_result") then
    open_dbg_context         = dump_multigrid_result
  else if (trim(context_name) .eq. "bsqvac_vac1") then
    open_dbg_context         = dump_bsqvac_vac1
  else if (trim(context_name) .eq. "phys_gc") then
    open_dbg_context         = dump_phys_gc
  else
    write(*,*) "unknown debug output context: '",trim(context_name),"'"
    stop
  end if

  ! check if debug out should be written at all
  should_write = iter2.le.max_dump
  open_dbg_context = open_dbg_context .and. should_write
  
  ! create output filename and open output file
  if (open_dbg_context) then
  
    ! TODO: debugging output into separate folder
  
    write(dump_filename, 995) trim(context_name), ns, iter2, trim(input_extension)
995 format(a,'_',i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(dump_filename)
  end if

end ! open_dbg_context

end module dbgout
