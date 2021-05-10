## DISCLAIMER

You are using a beta version of the PROGRAM VMEC, which is currently
under development by S. P. Hirshman at the Fusion Energy Division,
Oak Ridge National Laboratory.  Please report any problems or comments
to him.  As a beta version, this program is subject to change
and improvement without notice.

1. CODE SYNOPSIS

   THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
   SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
   FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
   REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
   ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
   LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
   VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
   DETERMINED BY MINIMIZING <M> = m\*\*2 S(m) , WHERE S(m) = Rm\*\*2 + Zm\*\*2 .
   AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
   NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
   A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
   USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
   VACUUM FIELD COMPONENTS BR, BPHI, BZ (see SUBROUTINE BECOIL)

   THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:

    B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) + iota(s) * grad(v) X grad(phiT)

   WHERE phiT is the toroidal flux (called phi in code) and
   u,v are the poloidal, toroidal angles, respectively.

2. ADDITIONAL CODES REQUIRED
   For the fixed boundary calculation, the user must provide the Fourier
   coefficients for the plasma boundary (the last surface outside of which
   the pressure gradient vanishes). For ALL but the simplest geometry, the
   SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
   code, can be used to produce the optimized VMEC Fourier representation for
   an arbritrary closed boundary (it need not be a 'star-like' DOmain, nor
   need it possess vertical, or 'stellarator', symmetry).

   For the free boundary calculation, the MAKEGRID code (available upon
   request) is needed to create a binary Green''s FUNCTION table for the
   vacuum magnetic field(s) and, IF data analysis is to be done, flux and
   field loops as well. The user provides a SUBROUTINE (BFIELD) which can be
   called at an arbitrary spatial location and which should RETURN the three
   cylindrical components of the vacuum field at that point. (Similary,
   locations of diagnostic flux loops, Rogowski coils, etc. are required IF
   equilibrium reconstruction is to be done.)

   Plotting is handled by a stand-alone package, PROUT.NCARG (written by
   R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
   file, WOUT.EXT, WHERE 'EXT' is the command-line extension of the INPUT file.

3. UNIX SCRIPT SETUP PARAMETERS
   The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
   the C-precompiler to produce both the machine-specific Fortran source and a
   make-file specific to ANY one of the following platforms:

   IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
   WINDOWS-NT, DEC-VMS

   Additional platforms are easy to add to the existing script as required.

4. FORTRAN PARAMETER STATEMENTS set by user
   In the Fortran-90 version of VMEC these PARAMETER statements have
   been replaced by dynamic memory allocation. So the user should set the
   run-time parameters ns (through ns_array), mpol, ntor in the NAMELIST INDATA.

### Added features since last edition
1. Implemented preconditioning algorithm for R,Z
2. The physical (unpreconditioned) residuals are used
   to determine the level of convergence
3. The original (MOMCON) scaling of lambda is used, i.e.,
   Bsupu = phip*(iota - lamda[sub]v)/SQRT(g). This is needed to
   maintain consistency with the time-stepper for arbitrary PHIP.

WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).

### Version History

 * 8.00
     1. added lforbal logical to fbal module to control whether to compute the flux-averaged
        force balance equation printed in the threed1 file. This requires a modification of
        the m=1,n=0 forces for R,Z in tomnsps subroutine. This works well, generally, and
        yields an improved <EQUIF> in threed1 file. However, there have been some cases where
        this non-variational departure fails to converge.
     2. added "bias" to iotas in bcovar when determining preconditioner for very small iota
        values. Seems to need this for improving convergence of preconditioned current-hole cases.
        Eliminated in v8.20.
     3. added "totzsp(s,a)_hess" to recompute r,z,l rapidly for only those modes that are jogged during
        Hessian calculation. NOTE: need to do this for lasym=true case, as well, eventually
 * 8.20  (August, 2004)
     1. removed 2-pt tie at origin for m=1 modes, to preserve tri-diagonal structure of Hessian.
        This is needed for preconditioning, which assumes block-tridi structure of equations
     2. fixed problem with free-boundary preconditioner, namely, ctor can not be extrapolated
        at edge when computing preconditioner, because this breaks tri-diagonal structure
     3. added new variables to input file to control preconditioning:
        1. PRECON_TYPE: = 'default', default tri-di (block size = 1)
                        = 'cg',      block tri-di, conjugate-gradient time-stepper
                        = 'gmres',   "          ", gmres time-stepper
                        = 'tfqmr',   "          ", transpose free qmr
        2. PREC2D_THRESHOLD: value of (unpreconditioned) forces at which block (2D) preconditioner
                             is turned on (=0 block preconditioner never turned on); recommended
                             (default) value ~ 1.E-10, or smaller, if convergence is poor
        3. LFORBAL: logical variable (default = .true.); when true, the force balance
                    used in the threed1 file is used to evolve the m=1 R,Z components. This
                    is a non-variational departure from the force equations for these modes,
                    but generally does not have an unfavorable impact on convergence.
     4. added new internal variable, ICTRL_PREC2D, to precon2d module. Replaces previous lprec2d
        and iequi>1 variables.
     5. removed lsweep_fast option from module precon2d. This slows the computation of the Hessian
        by about 2/3, but is more accurate (includes pdamp, liota, lforbal correctly)
     6. removed lflam logicals from bcovar and tomnsps, since now we must compute dFlam/dR,Z by
        jogging
     7. removed Compute_Hess_Flam_RZ from lamblks; this is now computed via jogging
        (also removed Get_dFlam_dRZ, FFT2Hessian, Forbal_avg, GetGsqrtVar supporting routines)
     8. removed internal liota logic, used to push iota profile rather than solving for it. Had
        needed this for symmetric Hessian (lsweep_fast=true option), but no longer required. Also,
        it was not implemented entirely correctly for lforbal=true case
     9. for lasym m=1 constraint rsc = zcc, changed xc array so that R+ = .5*(rsc + zcc) is stored at
        xc(rsc,m=1) and R- = .5*(rsc - zcc) is stored at xc(zcc,m=1). In residue, gcz(R-) == gcz(zcc)
        is zeroed by "evolving" gcr(zcc) = azd*[xc(zcc)-xcint], and gcr(rsc) => .5*[gcr(rsc) + gcz(zcc)]
        is evolved. In totzspa, the real rsc and zcc are computed from the internal representations
        (check convert call, too) by calling a new routine convert_asym (also called from wrout before
        writing out xc info). In readin, the original R+,R- are stored, so that for external "jogs",
        there will be no change in forces. All these changes are needed to obtain an invertible Hessian.
     10. added m=1 constraint for 3D case (similar to (i)), rss(n) = zcs(n), for n != 0. Imposed this
        on forces by adding routine constrain_m1 in residue. Added convert_sym routine to totzsp to convert
        from internal xc representation TO internal one.
     11. Decreased exponent on pdamp factor r2 (in bcovar) from 2 to 1, to give better conditioning
        especially for current hole cases
     12. Eliminated iotas bias for determining preconditioner, previously added in v8.00 for stabilizing
        current hole cases (not needed with corrected preconditioner)
 * 8.30  (October, 2004)
     1. Implemented flags for "reverse-communication" mode of vmec
 * 8.40
     1. Converted the m=1 constraints for 3D and asym back to old way; did not always
        converge well with the new constraints introduced in 8.20 (i-j)
 * 8.45  (December, 2005)
     1. Added the lconm1 logical. If = True, new constraint; if = False, old m=1 constraint used
     2. Added "perturbation" computation for lasym=TRUE case (totzspa_hess)
 * 8.46  (June, 2009)
     1. Added LRFP logical to allow easy switching of profiles between Stellarator/tokamak (PHIP=1, LRFP=F)
        and RFP (CHIP=1, LRFP=T). When LRFP=T, AI coefficients are expansion of q = 1/iota. Added lrfp to
        LIBSTELL/vmec_input module.
 * 8.47  (July, 2010)
     1. Rewrote magnetic field representation so that phip*lambda = new internal lambda. This greatly improves
        the conditioning of the lambda equations which otherwise become singular at the RFP reversal point
 * 8.48  (March 2012 - JDH)
     1. Accumulated small changes from SPH & JDH
     2. Modifications from J Geiger, March 2012
        - to be able to get additional main iterations if the force tolerance is
             not met. Parameter MAX_MAIN_ITERATIONS specifies how many main iteration
             cycles should be run.
        - to get a full output in the threed1-file if the force tolerance is not
             met. Specify the logical LFULL3D1OUT to true for this.
        - if vmec2000 is compiled with netcdf, you can still get the ascii-output
             if you specify the logical LWOUTTXT as true.
        - you get the output for diagno 1.0 and 1.5 if the logical LDIAGNO set true.
        - you get a rather old fort.8-output if you specify LOLDOUT as true.

        If none of these new variables is set, the behavior of vmec2000 is as
        expected from the version without the changes.
 * 8.49  (June, 2012)
     1. Fixed bug in bcovar when averaging half-grid covariant components onto full grid: did not
        zero components at (phantom) js=1 point, so edge force averages were incorrect
     2. Added lamscale factor to scale lambda in contravariant B-components. Modify wrout to maintain
        old-style lambda output
     3. Set lbsubs=F in jxbforce by default to capture current sheets
     4. Added lmove_axis INPUT logical (=T by default) so user can control whether or not the magnetic
        axis can be initially shifted to improve the initial force residuals. It is important NOT to move
        the helical axis for RFP equilibria requiring a helical seed (set l_moveaxis=F for this case!)

 * 8.50   (Jan, 2013)
     1. Improved scaling of lambda forces with respect to lamscale
     2. Fixed fnorm1 scaling (removed hs dependence)
     3. Added lgiveup logical (M. Drevlak/J. Geiger)

 * 8.51   (Sept, 2013)
     1. Restored hs-dependence of fnorm1 scaling (incorrectly) removed in 8.50

 * 8.52   (May, 2014)
     1. Fixed factor of 2 in integration factors tmult (WROUT) and dmult1 (JXBFORCE) for lasym=T,
        introduced when dnorm was changed in FIXARAY to include FULL theta range (1/23/14)
     2. In eqsolve (revert to v8.48), keep
           IF (liter_flag) CALL restart_iter(delt0)
        INSIDE IF (irst .eq. 2) TEST BLOCK
     3. Change default lforbal to FALSE in LIBSTELL, VMEC_INPUT (improves convergence,
        user can override in input file)

************************************************************************************
         
         HISTORY OF VMEC2000
 
************************************************************************************
 12/04:  Added INDATA namelist variables mfilter_bdy, nfilter_bdy to filter out high harmonics
         in the intial boundary representation for free boundary calculations. Default values (-1)
         DO NOT to filter ANY of the rbdy, zbdy harmonics. This is useful for very complex
         initial boundary shapes, for which it may be difficult to get the iterations started.
         Setting mfilter_bdy = M will exclude all harmonics in the initial boundary with m
         no. > M.
 11/04:  Added reverse communication loop in vmec.f, allowing runvmec to be called and controlled
         externally via an ictrl_array of flags.
 10/04:  Add GMRES, QMRS algorithms for stepping 2D, Preconditioner
 03/02:  Added LASYM to INDATA namelist to control running VMEC in symmetric (default) or
         asymmetric mode (stellarator symmetric...)
 03/01:  Introduced some damping on the free-boundary motion into SCALFOR in order to
         accelerate (and improve) the convergence of the free-bdy code.
 03/01:  Modified GUESS_AXIS subroutine to make better guess for axis. A grid search
         over the plasma cross section is made to find the axis value which maximizes
         the minimum value of the jacobian for each value of toroidal angle. This algorithm
         was implemented to improve convergence of VMEC2000 inside STELLOPT optimization code.
 
 03/00:  Begin modification of MHD forces to improve (post-processed) force balance at
         low aspect ratio. The first step is to divide out the factor of R in FR, FZ 
         forces. At low aspect ratio this produces mode-coupling which impacts the force.
         balance in the grad-s direction.
 01/00:  VMEC2000 is born!
         The latest modifications were made to improve the convergence with increasing
         mode numbers and finer radial meshes.
         Converted code to internal full mesh differencing for lambda, with a hybrid full/half
         differencing (blend in radius) for the lambda force differencing. The EXTRAP
         subroutine was eliminated, as well as the radial damping (in LAMCAL) on the
         lambda force. WROUT was modified so it continues to write out lambda on the
         half grid, consistent with previous versions of VMEC.
 
         Added (but temporarily disabled for backwards compatibility) a radial mesh zoning
         option to put the toroidal flux on a separate grid
         from the "s" grid. User can input aphi expansion coefficients for phi = phi(s)
         in the namelist: phi = aphi(1)*s + aphi(2)*s**2 ... The current, mass and iota
         expansions are now input as functions of the normalized toroidal flux (same as "s"
         in the old code). 
 
         Added an optional command line argument (look for "command line" below for further
         clarification) to specify the restart file name. The default, if lreset = F on the
         command line, is to look for the wout file with the input file name extension.
 
 08/99:  Fixed scaling of tcon with radial mesh (1000*hs**2),  to keep it invariant
 04/99:  Version 6.20: WROUT writes out signgs (VMEC now accepts pos jacobian)
         Additional diagnostic info is also output
 01/99:  Version 6.10: WROUT writes out pres, mass in pascals (1/mu0), 
         and currents in A(/m**2), NOT internal VMEC units.
 09/98:  Version 6.00. Start splitting out optimization routines.
         Decoupled BOOZER COORDINATE transformation from vmec/optimizer. Replaced with
         system call to xbooz
         Removed lspectrum_dump option. User can run xbooz code separately now.
         Removed loptim option. User now runs xstelopt code for optimization.
 07/98:  Added LOLDOUT logical to print out fort.8, fort.18 old-style output
         Added LEDGE_DUMP logical to print out fort.99 (edge-bfield) dump
         Added (1-cos(u))/2 weight to ripple in optimizer
 05/98:  Version 5.20 introduced as new marker. Added version ID in WOUT file to monitor changes.
         Fixed various subroutines to avoid 'segmentation fault' on DEC alphas,
         due to low stack space on those machines (eliminated large automatic
         vectors, replaced with dynamically allocated arrays). Users of DEC machines should
         check for adequate stack space using uname -a (or -s). If stack < 30000K, ask system
         administrator to increase it (although now VMEC will run with stack >= 2048K)
         Added Mercier criterion to optimization
         Added approximate external kink mode criterion to load_target
         Added Mercier condition calculation to jxbforce
         Fixed damping in lamcal: power = 0.5*m
         Added logical variable LMAC. If lmac = F (default), the
         mac-file (unit=nmac0) will be deleted (contains reconstruction data).
 04/98:  Added |B|mn-boozer spectral targets for to optimization chi-sq. (Version 5.10)
         Added LSPECTRUM_DUMP logical to indata file. If = T, will call load_target to dump
         boozer spectra, even for free/fixed boundary NOT being optimized.
         Writes to file bmn_spectrum.file-extension
         Fixed some 'bugs' associated with converting to a unique boundary representation
         used by the optimizer. Now the new input file written by the optimizer
         is consistent with the boundary representation
 03/98:  Improved boundary representation (unique_boundary, convert_boundary) for optimization
         code, as per Hirshman/Breslau representation.
 01/98:  Fixed default MGRID-FILE so user can specify ANY name, not necessarily
         prefixed by mgrid.    . (mgrid_file =      now acceptable in input file.)
 11/97:  Version 5.00: Fixed a number of F90-related bugs associated with
         the time-stepper and optimization loops. In particular, in EVOLVE, 
         when irst = 2, we now ALWAYS return before evolving XC, since GC 
         will NOT be the true force array coming out of FUNCT3D unless irst = 1.
         Also, the namelists have now been modularized (vmec_input, vmec_seq, optim_params modules).
         Added spres_ped: (0 < spres_ped <= 1) is a normed flux value. For
         s > spres_ped, the pressure profile is assumed flat. Default = 1 (no pedestal)
 08/97:  A separate sequence file, SEQ.EXT, is no longer supported. Now, to
         run a sequence of input files, add a namelist section &VSEQ to an
         input file, and the filenames (or extensions to input.ext..) will
         be sequentially executed.
         Added ntheta, nzeta to namelist (indata).
 08/97:  Added LPROF_OPT to optimizer namelist. When LPROF_OPT=T, then the AI (or AC)
         coefficients are allowed to vary. Thus, for ncurr=1, the current profile
         is varied, while for ncurr=0, the iota profile is varied (even though
         iota is matched in a chi-sq sense). When LPROF_OPT=F, AI (AC) coefficient array
         is fixed. LPROF_OPT replaces the older LCURPROF_OPT variable (which
         can still be read in, but is obsolete).
 07/97:  Beginning to phase out gcc/cc (c-precompilation). To this end, have introduced
         three new logical flags: LFREEB, LRECON, LOPTIM, which all default to
         F (fixed boundary, no reconstruction, no optimization) unless:
         (a) mgrid_file is specified in indata file, then LFREEB=T (unless LOPTIM=T is 
         specified);
         (b) itse or imse are nonzero in the indata file, then LRECON=T and
         LFREEB=T (unless LOPTIM=T or no mgrid_file exists)
         (c) if LOPTIM=T in indata, then LFREEB=LRECON=F regardless of their input
         values 
 07/97:  Established new makefile structure for VMEC. Now, the user can
         run vmec.lsqh script on a UNIX machine and only newer vmec files
         will be updated. To 'make' the vmec executable, type make debug (release)
         for a debug (release) version.
 04/97:  Began conversion of VMEC to F90 standard (Version 4.00). Making changes to
         argument passing conventions so as to be compatible with
         parallel (CRAY) processors. The code will no longer run under
         F77 with these changes. Converted all includes and common blocks
         to modules.
 07/96:  Added iresidue=3 condition to turn on FSQR(0,0) for fsq<FOPT_AXIS
         and turn off RADFOR-pfac time-variation
 07/96:  Moved pressum0 stuff into radfor routine
 07/96:  Fixed GETDIAM routine to use equilibrium pressure balance
         to compute diamagnetic flux correctly
 05/96:  Adjusted spatial damping parameter in FACLAM to improve convegence
 05/96:  Added constraint weighting-parameter, TCON0, to INDATA input file
 05/96:  WROUT modified: LMN Output on HALF-RADIAL mesh (same as internal VMEC mesh)
 05/96:  Added multigrid capability, NS_ARRAY and FTOL_ARRAY
 05/96:  Removed testing for MSE points changing sign in FIXRECON
 05/96:  Improved algorithm for moving axis to minimize CHISQ in AXISOPT
 04/96:  Upgraded fixed-boundary discrete p & iota profiles in profil1d
 04/96:  Remove Input_Update routine
 04/96:  Removed match to slope of MSE data in GETMSE routine
 04/96:  Added imatch_ip as imatch_edge=3 (for fixed boundary mode)
 03/96:  Made np (number field periods) a true parameter and added the
         variable nfper to vacuum include file.
 02/96:  Allow Pressure vs. s (s=normalized toroidal flux) input
         instead of Pressure vs. R (experimental input), when flag
         LPOFR = .FALSE. in namelist file. This is useful when TRANSP
         output, for example, is being used and the boundary shape is 
         uncertain.
 01/96:  ADDED OPTIONS FOR WINDOWS-NT COMPATABILITY AND DEBUGGING
         SPLIT FILES INTO SEPARATE .F MODULES
 11/95:  ASYMMETRIC THOUGHTS:
         We need to match Bpol = sqrt(BZ**2 + BR**2) in the
         MSE measurement, rather than just BZ!!! (Check with
         S.Batha,F.Levinton) 
 10/95:  Added Variable tension option for weighting tension
         on iota spline knots independently, i.e., namelist
         variables TENSI2, FPOLYI
 10/95:  Added NAMELIST variables ISNODES, IPNODES allowing
         user to pick number of iota,pressure knots
 10/95:  Accelerate the feedback in the imatch_phidege=2 loop
 09/95:  Merged INITSPLINE call into SETSPLINE routine
 06/95:  Converted input pressure coefficients (AM) to 
         be read in MKS units, NWT/M**2
 05/95:  Changed over to arbitrarily oriented B loops in GETBFLD
         Added namelist VSEQ for sequential running
         Added GETLIM, CAUCHY Subroutines to match to limiter position
 02/95:  Added Pres_Offset in Namelist for Thomson Data
 12/94:  Re-wrote GETTHOM, GETMSE subroutines to compute
         CHI-SQ directly from data (not SMOOTHED and SPLINED data).
         Added indexing arrays INDEX_I, INDEX_P in NEWPROFIL so
         calling order of GET... routines is now irrelevant.
         Looped around ALL GET... routines every IPEDSVD times.
         Computed extrapolated current in GETMSE correctly.
         Eliminated user options for fixed (in R) spline nodes
         in favor of fixed, uniformly distributed number of nodes
         in SQRT-S (or S) space.
         Rewrote CHISQ subroutine so all CHI-SQ's are computed
         from matrix elements, rather than from defining relations.
         Added VERSION NO. to THREED1 file output.
 10/94:  Improved pfac, phifac algorithms (a little)
 09/94:  Pass J-dot-B material through WOUT
      :  User should customize mgrid_defarea in subr readin for his site
 04/94:  LAMBDA DIFFERENCING ON HALF-MESH IMPLEMENTED
 03/94:  ELIMINATED RMNSS FOR PUSHING AXIS, REPLACED WITH RMNCC(JS=1)
         THIS WAS NECESSITATED BY DISCREPANCIES FOR SMALL ITSE...
 01/94:  IMPLEMENTED CHANGES FOR FIXED, FREE BOUNDARY
         FOR UP-DOWN NON-SYMMETRIC PLASMAS (symmetry_mode qualifier)
 11/93:  REPLACED WEIGHTS WITH STANDARD DEVIATIONS
 10/93:  IMPLEMENTED PHIEDGE MATCHING BASED ON WIDTH OF PRESSURE PROFILE DATA
 VERSION WHICH USES SQRT(S) MESH FOR SPLINES AND USES BOTH INNER AND
 OUTER EDGES FOR CURRENT MATCH (OPTIONAL RADIAL REDISTRIBUTION OF DATA)


## Input File Contents

               STANDARD INPUT DATA AND RECOMMENDED VALUES

  Plasma parameters (MKS units)
         ai:   expansion coefficients for iota (power series in s) used when ncurr=0
               Interpretation changes with piota_type
         am:   mass or pressure (gamma=0) expansion coefficients (series in s)
               in MKS units [NWT/M**2]
               Interpretation changes with pmass_type
         ac:   expansion coefficients for the normalized (pcurr(s=1) = 1)
               radial derivative of the flux-averaged toroidal current density
               (power series in s) used when ncurr=1
               Interpretation changes with pcurr_type
   ai_aux_s:   Auxiliary array for iota profile. Used for splines, s values
   ai_aux_f:   Auxiliary array for iota profile. Used for splines, function values
   am_aux_s:   Auxiliary array for mass profile. Used for splines, s values
   am_aux_f:   Auxiliary array for mass profile. Used for splines, function values
   ac_aux_s:   Auxiliary array for current profile. Used for splines, s values
   ac_aux_f:   Auxiliary array for current profile. Used for splines, function values
     curtor:   value of toroidal current [A]. Used if ncurr = 1 to specify
               current profile, or IF in data reconstruction mode.
    phiedge:   toroidal flux enclosed by plasma at edge (in Wb)
     extcur:   array of currents in each external current group. Used to
               multiply Green''s function for fields and loops read in from
               MGRID file. Should use real current units (A).
      gamma:   value of compressibility index (gamma=0 => pressure prescribed)
        nfp:   number of toroidal field periods ( =1 for Tokamak)
        rbc:   boundary coefficients of COS(m*theta-n*zeta) for R [m]
        zbs:   boundary coefficients of SIN(m*theta-n*zeta) for Z [m]
        rbs:   boundary coefficients of SIN(m*theta-n*zeta) for R [m]
        zbc:   boundary coefficients of COS(m*theta-n*zeta) for Z [m]


  Numerical and logical control parameters
      ncurr:   flux conserving (=0) or prescribed toroidal current (=1)
   ns_array:   array of radial mesh sizes to be used in multigrid sequence
   nvacskip:   number of iteration steps between accurate calculation of vacuum
               response; use fast interpolation scheme in between
 pres_scale:   factor used to scale pressure profile (default value = 1)
               useful so user can fix profile and change beta without having to change
               all AM coefficients separately
      tcon0:   weight factor for constraint force (=1 by DEFAULT)
      lasym:   =T, run in asymmetric mode; =F, run in stellarator symmetry mode
     lfreeb:   =T, run in free boundary mode if mgrid_file exists
    lforbal:   =T, use non-variational forces to ensure <EQUIF> = 0;
               =F, use variational form of forces, <EQUIF> ~ 0

  Convergence control parameters
 ftol_array:   array of value of residual(s) at which each multigrid
               iteration ends
niter_array:   array of number of iterations (used to terminate run) at
               each multigrid iteration
      nstep:   number of timesteps between printouts on screen
   nvacskip:   iterations skipped between full update of vacuum solution

  Preconditioner control parameters (added 8/30/04)
precon_type:   specifies type of 2D preconditioner to use:
               'default'            - diagonal in m,n, tri-diagonal in s
               'conjugate-gradient' - block tri-di, evolve using cg method
               'gmres'              - block tri-di, generalized minimal residual method
               'tfqmr'              - block tri-di, transpose-free quasi minimum residual
prec2d_threshold:
               value of preconditioned force residuals at which block (2d) tri-di
               solver is turned on, if requested via type_prec2d

  Character parameters
 mgrid_file:   full path for vacuum Green''s function data
 pcurr_type:   Specifies parameterization type of pcurr function
                 'power_series' - I'(s)=Sum[ ac(j) s ** j] - Default
                 'gauss_trunc'  - I'(s)=ac(0) (exp(-(s/ac(1)) ** 2) -
                                               exp(-(1/ac(1)) ** 2))
                  others - see function pcurr
 piota_type:   Specifies parameterization type of piota function
                 'power_series' - p(s)=Sum[ am(j) s ** j] - Default
                  others - see function piota
 pmass_type:   Specifies parameterization type of pmass function
                 'power_series' - p(s)=Sum[ am(j) s ** j] - Default
                 'gauss_trunc'  - p(s)=am(0) (exp(-(s/am(1)) ** 2) -
                                               exp(-(1/am(1)) ** 2))
                  others - see function pmass
