module dbgout
  use json
  implicit none

  logical :: skip_dbgout_collison = .false.

contains

!> check if any output is desired for the current iteration
!> check if the given context should be openend based on input file flags
!> check if
!> @param context_name a string describing the subroutine from which this function is called
!> @param repetition   a number to distinguish two calls to this with the same value of iter2
!> @param id           a number to replace iter2 in the output filename
!> @return true: debug output should be written and file is open; false otherwise
function open_dbg_context(context_name, repetition, id)
  use vmec_dim,   only: ns
  use vmec_main,  only: iter2
  use vmec_input
  implicit none

  character(len=*), intent(in)  :: context_name
  integer, intent(in), optional :: repetition
  integer, intent(in), optional :: id
  logical :: open_dbg_context

  character(len=255) :: dump_filename
  character(len=255) :: output_folder
  logical            :: should_write, file_exists
  integer            :: iter_value_to_use, i

  ! enable semi-pretty-printing JSON data
  json_pretty_print = .true.

  if (present(id)) then
    iter_value_to_use = id
  else
    iter_value_to_use = iter2
  end if

  ! check if debug out should be written at all
  should_write = .false.
  do i = 1, num_iter2_to_dump
    if (iter_value_to_use .eq. iter2_to_dump(i)) then
      should_write = .true.
    end if
  end do ! num_iter2_to_dump

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
  else if (trim(context_name) .eq. "totzsp_input") then
    open_dbg_context         = dump_totzsp_input
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
  else if (trim(context_name) .eq. "phys_gc") then
    open_dbg_context         = dump_phys_gc
  else if (trim(context_name) .eq. "multigrid_result") then
    open_dbg_context         = dump_multigrid_result

    ! multigrid_result needs to be written once at end of many iterations,
    ! so the usual should_write logic needs to be broken here
    should_write = .true.

  else if (trim(context_name) .eq. "printout") then
    open_dbg_context         = dump_printout

  ! fileout
  else if (trim(context_name) .eq. "bcovar_fileout") then
    open_dbg_context         = dump_bcovar_fileout
  else if (trim(context_name) .eq. "bss") then
    open_dbg_context         = dump_bss
  else if (trim(context_name) .eq. "jxbforce_bsub_lowpass") then
    open_dbg_context         = dump_jxbforce_bsub_lowpass
  else if (trim(context_name) .eq. "jxbout") then
    open_dbg_context         = dump_jxbout
  else if (trim(context_name) .eq. "mercier") then
    open_dbg_context         = dump_mercier
  else if (trim(context_name) .eq. "threed1_firstTable") then
    open_dbg_context         = dump_threed1_firstTable
  else if (trim(context_name) .eq. "threed1_geomag") then
    open_dbg_context         = dump_threed1_geomag
  else if (trim(context_name) .eq. "threed1_volquant") then
    open_dbg_context         = dump_threed1_volquant
  else if (trim(context_name) .eq. "threed1_axis") then
    open_dbg_context         = dump_threed1_axis
  else if (trim(context_name) .eq. "threed1_beta") then
    open_dbg_context         = dump_threed1_beta
  else if (trim(context_name) .eq. "threed1_shafrint") then
    open_dbg_context         = dump_threed1_shafrint
  else if (trim(context_name) .eq. "freeb_data") then
    open_dbg_context         = dump_freeb_data


  ! NESTOR vac1
  else if (trim(context_name) .eq. "vac1n_vacuum") then
    open_dbg_context         = dump_vac1n_vacuum
  else if (trim(context_name) .eq. "vac1n_precal") then
    open_dbg_context         = dump_vac1n_precal
  else if (trim(context_name) .eq. "vac1n_surface") then
    open_dbg_context         = dump_vac1n_surface
  else if (trim(context_name) .eq. "vac1n_bextern") then
    open_dbg_context         = dump_vac1n_bextern
  else if (trim(context_name) .eq. "vac1n_analyt") then
    open_dbg_context         = dump_vac1n_analyt
  else if (trim(context_name) .eq. "vac1n_greenf") then
    open_dbg_context         = dump_vac1n_greenf
  else if (trim(context_name) .eq. "vac1n_fourp") then
    open_dbg_context         = dump_vac1n_fourp
  else if (trim(context_name) .eq. "vac1n_fouri") then
    open_dbg_context         = dump_vac1n_fouri
  else if (trim(context_name) .eq. "vac1n_solver") then
    open_dbg_context         = dump_vac1n_solver
  else if (trim(context_name) .eq. "vac1n_bsqvac") then
    open_dbg_context         = dump_vac1n_bsqvac

  ! NESTOR vac2
  else if (trim(context_name) .eq. "vac2_vacuum") then
    open_dbg_context         = dump_vac2_vacuum
  else if (trim(context_name) .eq. "vac2_precal") then
    open_dbg_context         = dump_vac2_precal
  else if (trim(context_name) .eq. "vac2_surface") then
    open_dbg_context         = dump_vac2_surface
  else if (trim(context_name) .eq. "vac2_bexmat") then
    open_dbg_context         = dump_vac2_bexmat
  else if (trim(context_name) .eq. "vac2_matrix") then
    open_dbg_context         = dump_vac2_matrix
  else if (trim(context_name) .eq. "vac2_foumat_unreg") then
    open_dbg_context         = dump_vac2_foumat_unreg
  else if (trim(context_name) .eq. "vac2_analin") then
    open_dbg_context         = dump_vac2_analin
  else if (trim(context_name) .eq. "vac2_analyt") then
    open_dbg_context         = dump_vac2_analyt
  else if (trim(context_name) .eq. "vac2_foumat") then
    open_dbg_context         = dump_vac2_foumat
  else if (trim(context_name) .eq. "vac2_linsys") then
    open_dbg_context         = dump_vac2_linsys
  else if (trim(context_name) .eq. "vac2_linslv") then
    open_dbg_context         = dump_vac2_linslv
  else if (trim(context_name) .eq. "vac2_bsqvac") then
    open_dbg_context         = dump_vac2_bsqvac


  ! default
  else
    write(*,*) "unknown debug output context: '",trim(context_name),"'"
    stop
  end if

  open_dbg_context = open_dbg_context .and. should_write

  ! create output filename and open output file
  if (open_dbg_context) then

    ! debugging output into separate folder "input_extension"
    output_folder = trim(input_extension) // "/" // trim(context_name)
    CALL system("mkdir -p "//trim(output_folder)) ! NOTE: This only works on Linux/Unix !!!

    if (present(id)) then
      if (present(repetition)) then
        write(dump_filename, 998) trim(output_folder), &
                                  trim(context_name),  &
                                  ns, id, repetition, &
                                  trim(input_extension)
      else
        write(dump_filename, 999) trim(output_folder), &
                                  trim(context_name),  &
                                  ns, id, trim(input_extension)
      end if
    else
      if (present(repetition)) then
        write(dump_filename, 998) trim(output_folder), &
                                  trim(context_name),  &
                                  ns, iter2, repetition, &
                                  trim(input_extension)
      else ! default: ns, iter2 for filename
        write(dump_filename, 999) trim(output_folder), &
                                  trim(context_name),  &
                                  ns, iter2, trim(input_extension)
      end if
    end if
998   format(a,'/',a,'_',i5.5,'_',i6.6,'_',i2.2,'.',a,'.json')
999   format(a,'/',a,'_',i5.5,'_',i6.6,'_01.',a,'.json')

    ! check if file already exists (and stop in that case)
    inquire(file=trim(dump_filename), exist=file_exists)
    if (file_exists) then
      if (skip_dbgout_collison) then
        ! Temporary hack to skip overwriting a file
        ! without halting the algorithm.
        return
      end if
      stop "debug output file already exists: '"//trim(dump_filename)//"'"
    end if

    call open_dbg_out(dump_filename)
  end if

end ! open_dbg_context

end module dbgout
