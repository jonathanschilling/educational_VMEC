!-----------------------------------------------------------------------
!     Module:        read_spec_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          03/26/2012
!     Description:   This subroutine reads the SPEC HDF5 file.  In order
!                    to ease integration into existing routines the
!                    module initializes the varialbes in the
!                    read_wout_mod VMEC output interface module.
!-----------------------------------------------------------------------
      MODULE read_spec_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE read_wout_mod, ONLY: nfp, ns, mpol,ntor, mnmax, rmax_surf, &
                               rmin_surf, zmax_surf, Rmajor, Itor, &
                               betatot, &
                               rmnc, rmns, zmnc, zmns, lmnc, lmns, &
                               bsupumnc, bsupumns, bsupvmnc, bsupvmns, &
                               iotaf, presf, phi, xm, xn, lasym, &
                               lthreed, input_extension, version_
!-----------------------------------------------------------------------
!     Local Variables
!          ier            Error Flag
!          iunit          File ID Number
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PRIVATE :: ier, i, j, k
!-----------------------------------------------------------------------
!     Subroutines
!          read_spec_file     Read the SPEC output file.
!-----------------------------------------------------------------------
      
      INTERFACE read_spec_file
          MODULE PROCEDURE read_spec_hdf5
      END INTERFACE
      
      CONTAINS
      
      SUBROUTINE read_spec_hdf5(id_string)
      CHARACTER(LEN=*), INTENT(in) :: id_string
      INTEGER :: ier, dex
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
      ! Check for extension or full name
      input_extension = ''
      STOP 'No HDF5 Support, cannot read SPEC output'
!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------    
      END SUBROUTINE read_spec_hdf5
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE read_spec_mod
