module dbgout
  use json

  implicit none

  character(len=255) :: dump_filename

  ! maximum number of iterations for which to dump data
  integer, parameter :: max_dump = 10

  logical, parameter :: dump_add_fluxes          = .false.
  logical, parameter :: dump_metric              = .false.
  logical, parameter :: dump_volume              = .false.
  logical, parameter :: dump_bcontrav            = .false.
  logical, parameter :: dump_bcov                = .false.
  logical, parameter :: dump_lambda_forces       = .false.
  logical, parameter :: dump_bcov_full           = .false.
  logical, parameter :: dump_precondn            = .false.
  logical, parameter :: dump_forceNorms_tcon     = .false.
  logical, parameter :: dump_lulv_comb           = .false.
  logical, parameter :: dump_calc_fbal           = .false.
  logical, parameter :: dump_evolve              = .false.
  logical, parameter :: dump_fixaray             = .false.
  logical, parameter :: dump_spectral_constraint = .false.
  logical, parameter :: dump_forces              = .false.
  logical, parameter :: dump_geometry            = .false.
  logical, parameter :: dump_constraint_force    = .false.
  logical, parameter :: dump_guess_axis          = .false.
  logical, parameter :: dump_interp              = .false.
  logical, parameter :: dump_jacobian            = .false.
  logical, parameter :: dump_lamcal              = .false.
  logical, parameter :: dump_profil1d            = .false.
  logical, parameter :: dump_profil3d            = .false.
  logical, parameter :: dump_readin_boundary     = .false.
  logical, parameter :: dump_physical_gc         = .false.
  logical, parameter :: dump_fsq                 = .false.
  logical, parameter :: dump_scale_m1            = .false.
  logical, parameter :: dump_scalfor_out         = .false.
  logical, parameter :: dump_fsq1                = .false.
  logical, parameter :: dump_scalfor             = .false.
  logical, parameter :: dump_symforce            = .false.
  logical, parameter :: dump_tomnsps             = .false.
  logical, parameter :: dump_tomnspa             = .false.
  logical, parameter :: dump_multigrid_result    = .false.
  logical, parameter :: dump_bsqvac              = .false.
  logical, parameter :: dump_rbsq                = .false.

end module dbgout
