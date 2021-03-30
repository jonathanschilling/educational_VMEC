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
