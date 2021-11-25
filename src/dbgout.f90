module dbgout
  use json

  implicit none

  character(len=255) :: dump_filename

  ! maximum number of iterations for which to dump data
  integer :: max_dump                 = 10

  logical :: dump_add_fluxes          = .false.
  logical :: dump_metric              = .false.
  logical :: dump_volume              = .false.
  logical :: dump_bcontrav            = .false.
  logical :: dump_bcov                = .false.
  logical :: dump_lambda_forces       = .false.
  logical :: dump_bcov_full           = .false.
  logical :: dump_precondn            = .false.
  logical :: dump_forceNorms_tcon     = .false.
  logical :: dump_lulv_comb           = .false.
  logical :: dump_calc_fbal           = .false.
  logical :: dump_evolve              = .false.
  logical :: dump_fixaray             = .false.
  logical :: dump_spectral_constraint = .false.
  logical :: dump_forces              = .false.
  logical :: dump_geometry            = .false.
  logical :: dump_constraint_force    = .false.
  logical :: dump_guess_axis          = .false.
  logical :: dump_interp              = .false.
  logical :: dump_jacobian            = .false.
  logical :: dump_lamcal              = .false.
  logical :: dump_profil1d            = .false.
  logical :: dump_profil3d            = .false.
  logical :: dump_readin_boundary     = .false.
  logical :: dump_physical_gc         = .false.
  logical :: dump_fsq                 = .false.
  logical :: dump_scale_m1            = .false.
  logical :: dump_scalfor_out         = .false.
  logical :: dump_fsq1                = .false.
  logical :: dump_scalfor             = .false.
  logical :: dump_symforce            = .false.
  logical :: dump_tomnsps             = .false.
  logical :: dump_tomnspa             = .false.
  logical :: dump_multigrid_result    = .false.
  logical :: dump_bsqvac              = .false.
  logical :: dump_rbsq                = .false.

end module dbgout
