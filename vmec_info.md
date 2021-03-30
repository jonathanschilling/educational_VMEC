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
     3. added "totzsps,a_hess" to recompute r,z,l rapidly for only those modes that are jogged during
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

## Input File Contents

      LOCAL VARIABLES

      rbcc,rbss,rbcs,rbsc
               boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
      zbcc,zbss,zbcs,zbsc
               boundary Fourier coefficient arrays for Z

      XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC

      STACKING ORDER DEPENDS ON LASYM AND LTHREED. EACH COMPONENT XCC, XSS, XSC, XCS
      HAS SIZE = mns. (PHIFAC, MSE TAKE UP 1 INDEX EACH AT END OF ARRAY)

        LTHREED=F,      LTHREED=F,      LTHREED=T,      LTHREED=T
        LASYM=F         LASYM=T         LASYM=F         LASYM=T

         rmncc           rmncc           rmncc           rmncc
         zmnsc           rmnsc           rmnss           rmnss
         lmnsc           zmnsc           zmnsc           rmnsc
                         zmncc           zmncs           rmncs
                         lmnsc           lmnsc           zmnsc
                         lmncc           lmncs           zmncs
                                                         zmncc
                                                         zmnss
                                                         lmnsc
                                                         lmncs
                                                         lmncc
                                                         lmnss


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
precon_type:   specifies type of 2D preconditioner to use ('default', diagonal in m,n,
               tri-diagonal in s; 'conjugate-gradient', block tri-di, evolve using
               cg method; 'gmres', block tri-di, generalized minimal residual method;
               'tfqmr', block tri-di, transpose-free quasi minimum residual
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

  Equilibrium reconstruction parameters
     phifac:   factor scaling toroidal flux to match apres or limiter
  datastark:   pitch angle data from stark measurement
   datathom:   pressure data from Thompson, CHEERS (Pa)
    imatch_         = 1 (default),match value of PHIEDGE in input file
    phiedge:   = 0, USE pressure profile width to determine PHIEDGE
               = 2, USE LIMPOS data (in mgrid file) to find PHIEDGE
               = 3, USE Ip to find PHIEDGE (fixed-boundary only)
       imse:   number of Motional Stark effect data points
               >0, USE mse data to find iota; <=0, fixed iota profile ai
       itse:   number of pressure profile data points
               = 0, no thompson scattering data to READ
    isnodes:   number of iota spline points (computed internally unless specified explicitly)
    ipnodes:   number of pressure spline points (computed internally unless specified explicitly)
      lpofr:   LOGICAL variable. =.true. IF pressure data are
               prescribed in REAL space. =.false. IF data in flux space.
     pknots:   array of pressure knot values in SQRT(s) space
     sknots:   array of iota knot values in SQRT(s) space
      tensp:   spline tension for pressure profile

      tensi:   spline tension for iota
     tensi2:   vbl spline tension for iota
     fpolyi:   vbl spline tension form factor (note: IF tensi!=tensi2
              THEN tension(i-th point) = tensi+(tensi2-tensi)*(i/n-1))**fpolyi
              - - - - - - - - - - - - - - - - - -
   mseangle_   uniform EXPerimental offset of MSE data
    offset:    (calibration offset) ... PLUS ...
   mseangle_   multiplier on mseprof offset array
    offsetM:   (calibration offset)
    mseprof:   offset array from NAMELIST MSEPROFIL
               so that the total offset on the i-th MSE data point is
               taken to be
               = mseangle_offset+mseangle_offsetM*mseprof(i)
              - - - - - - - - - - - - - - - - - -
pres_offset:   uniform arbitrary  radial offset of pressure data
    presfac:   number by which Thomson scattering data is scaled
               to get actual pressure
    phidiam:   diamagnetic toroidal flux (Wb)
     dsiobt:   measured flux loop signals corresponding to the
               combination of signals in iconnect array
    indxflx:   array giving INDEX of flux measurement in iconnect array
   indxbfld:   array giving INDEX of bfield measurement used in matching
       nobd:   number of connected flux loop measurements
     nobser:   number of individual flux loop positions
     nbsets:   number of B-coil sets defined in mgrid file
 nbcoils(n):   number of bfield coils in each set defined in mgrid file
   nbcoilsn:   total number of bfield coils defined in mgrid file
   bbc(m,n):   measured magnetic field at rbcoil(m,n),zbcoil(m,n) at
               the orientation br*COS(abcoil) + bz*SIN(abcoil)
rbcoil(m,n):   R position of the m-th coil in the n-th set from mgrid file
zbcoil(m,n):   Z position of the m-th coil in the n-th set from mgrid file
abcoil(m,n):   orientation (surface normal wrt R axis; in radians)
               of the m-th coil in the n-th set from mgrid file.
      nflxs:   number of flux loop measurements used in matching
   nbfld(n):   number of selected EXTERNAL bfield measurements in set n from nml file
     nbfldn:   total number of EXTERNAL bfield measurements used in matching
              - - - - - - - - - - - - - - - - - -
            NOTE: FOR STANDARD DEVIATIONS (sigma''s) < 0, INTERPRET
            AS PERCENT OF RESPECTIVE MEASUREMENT
 sigma_thom:   standard deviation (Pa) for pressure profile data
sigma_stark:   standard deviation (degrees) in MSE data
 sigma_flux:   standard deviaton (Wb) for EXTERNAL poloidal flux data
    sigma_b:   standard deviation (T) for EXTERNAL magnetic field data
sigma_current:  standard deviation (A) in toroidal current
sigma_delphid:  standard deviation (Wb) for diamagnetic match


      THE (ABSOLUTE) CHI-SQ ERROR IS DEFINED AS FOLLOWS:

         2
      CHI      =     SUM [ EQ(K,IOTA,PRESSURE)  -  DATA(K) ] ** 2
                    (K) -----------------------------------
                                  SIGMA(K)**2

      HERE, SIGMA IS THE STANDARD DEVIATION OF THE MEASURED DATA, AND
      EQ(IOTA,PRESSURE) IS THE EQUILIBRIUM EXPRESSION FOR THE DATA TO BE
      MATCHED:

      EQ(I)   =    SUM [ W(I,J)*X(J) ]
                  (J)

      WHERE W(I,J) ARE THE (LINEAR) MATRIX ELEMENTS AND X(J) REPRESENT
      THE KNOT VALUES OF IOTA (AND/OR PRESSURE). THE RESULTING LEAST-SQUARES
      MATRIX ELEMENTS AND DATA ARRAY CAN BE EXPRESSED AS FOLLOWS:

      ALSQ(I,J) = SUM [ W(K,I) * W(K,J) / SIGMA(K) ** 2]
                  (K)

      BLSQ(I)   = SUM [ W(K,I) * DATA(K)/ SIGMA(K) ** 2]
                  (K)

      THEREFORE, INTERNALLY IT IS CONVENIENT TO WORK WITH THE 'SCALED'
      W'(K,I) = W(K,I)/SIGMA(K) AND DATA'(K) = DATA(K)/SIGMA(K)

      ****!   I - M - P - O - R - T - A - N - T     N - O - T - E   *****

      THE INPUT DATA FILE WILL ACCEPT BOTH POSITIVE AND NEGATIVE
      SIGMAS, WHICH IT INTERPRETS DIFFERENTLY. FOR SIGMA > 0, IT
      TAKES SIGMA TO BE THE STANDARD DEVIATION FOR THAT MEASUREMENT
      AS DESCRIBED ABOVE. FOR SIGMA < 0, SIGMA IS INTERPRETED AS
      THE FRACTION OF THE MEASURED DATA NEEDED TO COMPUTE THE ABSOLUTE
      SIGMA, I.E., (-SIGMA * DATA) = ACTUAL SIGMA USED IN CODE.
