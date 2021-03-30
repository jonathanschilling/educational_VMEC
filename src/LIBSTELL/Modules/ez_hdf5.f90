!-----------------------------------------------------------------------
!     Module:        ez_hdf5
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          12/5/2012
!     Description:   This module helps simplfy the interface to hdf5
!                    files.  The HDF5 module must be loaded for your
!                    compiler. The following compiler commands should be
!                    added to your makefiles to properly compile this
!                    file:
!
!                    -L$(HDF5_HOME)/lib -I$(HDF5_HOME)/include 
!                    -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 
!                    -lpthread -lz -lm
!
!                    At the end of this file there is a commented
!                    section which contains an example showing how to
!                    call these routines.  
!-----------------------------------------------------------------------
      MODULE ez_hdf5
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE ez_hdf5
!
!
!
!-----  B E G I N  S A M P L E  C O D E  -------------------------------
!      SUBROUTINE WRITE_HDF5_TEST
!      USE ez_hdf5
!      IMPLICIT NONE
!      INTEGER :: ier
!      REAL    :: test_arr(4,6)
!-----------------------------------------------------------------------
!     Begin Subroutine
!-----------------------------------------------------------------------
!      test_arr(:,:) = 1.56
!      !-----  OPEN HDF5 FILE
!      PRINT *,'WRITING TO FILE: ','test.h5'
!      CALL open_hdf5('test.h5',fid,ier)
!      !-----  WRITE SCALARS
!      CALL write_var_hdf5(fid,'int_var',ier,INTVAR=314,ATT='An Integer Variable',ATT_NAME='Description')
!      CALL write_var_hdf5(fid,'flt_var',ier,FLTVAR=3.14159,ATT='Approximation of Pi')
!      CALL write_var_hdf5(fid,'boo_var',ier,BOOVAR=.TRUE.,ATT='A Boolean Variable')
!      !-----  WRITE VECTORS
!      CALL write_var_hdf5(fid,'int_arr',5,ier,INTVAR=/0,1,2,3,4/)
!      !-----  WRITE ARRAYS
!      CALL write_var_hdf5(fid,'rmnc',4,6,ier,FLTVAR=test_arr,ATT='4x6 Test Array')
!      !-----  CLOSE HDF5 FILE
!      CALL close_hdf5(fid,ier)
!      END SUBROUTINE WRITE_HDF5_TEST
!-----------------------------------------------------------------------
