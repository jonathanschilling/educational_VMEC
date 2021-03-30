!SAME as adas_mog.f90 but have to be sure that all tables existin  $ADASDIR/tables...
!---------------------------------------------------------------------
! adas_mog.f90: by M.Gorelenkova                                      
!               created 07/2009                                       
!               based on Michael Kraus routines                       
!                                                                     
!  Module which computes reaction rate coefficient tables once and stores these as files;
!  once a table is computed, it is not recomputed, but simply read in
!  for interpolation.  
!  All energies & temperatures in KeV/amu
!  All rate coefficients <sigma*v> in m**3/sec.
!
! subroutine adas_btsigv -- calculates beam-target maxwellian average <sig*v> (Eb/Ab,Ti/Ai)
!                           for neutralizing CX and II reactions, available for H and He. 

!
! subroutine adas_adas_sigvte_ioniz  --  calculates ionization of neutrals by electron impact
!                                        rate coefficients <sigma*v> averaged over Maxwellian electron
!                                        distribution characterized by temperature Te (KeV),
!                                        available only for H. 
!
! subroutine adas_zstop_sigvz -- calculates sig*v tables of neutral stopping   
!                                on fully stripped light impurities (vs. Erel only); 
!                                as a sum of sig*v for ionization (includes non-neutralazing CX). 
!                                Data for H and He neutral atoms are available.
!
! subroutine adas_sigv -- calculates sig*v(Eb/Ab) for neutralizing charge-exchange reaction 
!                         and impact ionization which includes  non-neutralizing 
!                         charge-exchange reaction. Data for H and He neutral atoms are available.
!
! subroutine adas_sigv --  creates tables for  CX and ionization cross section sig(Eb/Ab). 
!                          Data for H and He neutral atoms are available. 
!
! subroutine adas_bms  -- computes beam stopping rate coefficient on impurities for H - beam.
!                         Routine uses  EQUIVALENT of  Electron Density [cm**-3]
!                         N_el*SUM(Z_imp^2* N_imp)/SUM(Z_imp * N_imp)/Z_imp, where
!                         N_el -- electron density and impurity density 
!                         N_imp -- impurity density 
!                         for details see [1],p.794
!
! All routines try to open pre-computed tables for an asking reaction in $ADASDIR/tables/...
! directory. If these tables had been created early then data will be read and interpolated, if
! not then table will be computed and wrote to the $ADASDIR/tables/... first.
!
! References:                                                         
!                                                                     
!    1. H. Anderson et al. (2000)                                     
!       "Neutral beam stopping and emission in fusion plasmas I"      
!       Plasma Physics and Controlled Fusion Vol. 42, pp 781-806      
!                                                                     
!    2. H. P. Summers (2004)                                          
!       "The ADAS User Manual v2.6"                                   
!       http://www.adas.ac.uk/manual.php                              
!                                                                     
!---------------------------------------------------------------------
Module adas_mod_simpl

end module adas_mod_simpl
