!
!  This program developed in FORTRAN is a Finite Element solver for linear-
!  static analyses as well as linearized stability analyses. It is inherently
!  coupled to the open-source panel method APAME for providing fluid-structure-
!  interaction capabilites.
!    
!  Copyright (C) 2024 TUD Dresden University of Technology
! 
!  This file is part of FiPPS².
!
!  FiPPS² is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  FiPPS² is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
! =================================================================================================
!
!> @brief
!> Subroutine estimates the number of nonzero-enries in global geometric
!> stiffness matrix using special matrix structure of element geometric
!> stiffness matrix
!> This number is needed for sufficient memory allocation in PETSc.
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter, 06.05.2010
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 22.06.2010
!
!> $Id: sol2_extract_eigenvalues.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol2_extract_eigenvalues(eps, eigvals, eigenvectors, nconv, dim_dof, rank)
!
! use
  use konstanten
!
! =================================================================================================
!
! Include
!
#include "slepc/finclude/slepceps.h"
  use slepceps
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  integer                                         :: ii
  integer                                         :: dim_dof
  
  double precision, dimension(nconv), intent(out) :: eigvals
  double precision, dimension(dim_dof,nconv), intent(out) :: eigenvectors
  
  
  Vec                                             :: eigvec       ! Eigenvektor als Sparse-Vektor
  Vec                                             :: eigvec_seq
  EPS                                             :: eps
  PetscScalar                                     :: kr, ki
  PetscInt                                        :: nconv
  PetscErrorCode                                  :: ierr
  PetscInt                                        :: pii
  PetscScalar, pointer                            :: xx_v(:)
  PetscMPIInt                                     :: rank
  PetscReal                                       :: error
  VecScatter                                      :: scatter
  
  integer                                         :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Display solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!     ** Get number of converged eigenpairs
      if (textoutput .eq. .true. .and. rank .eq. 0) write(*,150) nconv
 150  format (' Number of converged eigenpairs:',I2/)
!     ** Display eigenvalues and relative errors
      if (nconv.gt.0) then
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '         k          ||Ax-kx||/||kx||'
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) ' ----------------- ------------------'
        call VecCreate(PETSC_COMM_WORLD,eigvec,ierr); CHKERRQ(ierr)
        call VecSetSizes(eigvec,PETSC_DECIDE,dim_dof,ierr); CHKERRQ(ierr)
        call VecSetFromOptions(eigvec,ierr); CHKERRQ(ierr)
 
        do pii=0, nconv-1
!          ** Get converged eigenpairs: i-th eigenvalue is stored in kr 
!          ** (real part) and ki (imaginary part)
           call EPSGetEigenpair(eps,pii,kr,ki,eigvec,PETSC_NULL_VEC,ierr); CHKERRQ(ierr)
           
           call VecCreate(PETSC_COMM_WORLD,eigvec_seq,ierr); CHKERRQ(ierr)
           
           call VecScatterCreateToZero(eigvec, scatter, eigvec_seq, ierr); CHKERRQ(ierr)
   
           call VecScatterBegin(scatter, eigvec, eigvec_seq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
           call VecScatterEnd  (scatter, eigvec, eigvec_seq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    
           call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
           
           if (rank .eq. 0) then

             call VecGetArrayF90(eigvec_seq,xx_v,ierr); CHKERRQ(ierr)
             ! Access first local entry in vector with
             do ii = 1, dim_dof
               eigenvectors(ii,pii+1) = xx_v(ii)
             end do

             call VecRestoreArrayF90(eigvec_seq,xx_v,ierr); CHKERRQ(ierr)
 
             eigvals(pii+1) = PetscRealPart(kr)
           end if
           
           call VecDestroy(eigvec_seq, ierr); CHKERRQ(ierr)
           
           ! ** Compute the relative error associated to each eigenpair
           if (textoutput .eq. .true.) then 
             call EPSComputeError(eps,pii,EPS_ERROR_RELATIVE,error,ierr); CHKERRQ(ierr)
             if (rank .eq. 0) write(*,160) PetscRealPart(kr), error
 160         format (1P,'   ',E12.4,'     ',E12.4)
           end if
        enddo
        
        if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*)
        call VecDestroy(eigvec, ierr); CHKERRQ(ierr)
        
      endif
   
      ! calculate reciprocal eigenvalue and switch algebraic sign
      ! reciprocal since SLEPc calculate the eigenvalue problem Kg*x=mu*K*x
      ! switch in algebraic sign since membrane forces are treated as N=lambda*N0
      ! instead of N=-lambda*N0 in Routine for geometric stiffness matrix
      
      if (rank .eq. 0) then
        do ii=1,nconv

          eigvals(ii)=-1.D0/eigvals(ii)

        end do
      end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol2_extract_eigenvalues'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine sol2_extract_eigenvalues
