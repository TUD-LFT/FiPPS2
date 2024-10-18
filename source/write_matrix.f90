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
!> Subroutine outputs a PETSc-Matrix
!
!> @details
!
!> @author 
!
!> $Id: write_matrix.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine write_matrix(K, name)
  
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscmat.h"
  use petscmat
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  
  Mat            :: K
  PetscViewer    :: viewer
  PetscErrorCode :: ierr
  
  character(*) :: name
  
  integer               :: err_code=0
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
  !call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER,PETSC_NULL_CHARACTER,0,0,1800,1150,viewer,ierr); CHKERRQ(ierr)
  !call PetscViewerBinaryOpen(PETSC_COMM_WORLD, TRIM(name), FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerASCIIOpen(PETSC_COMM_WORLD, TRIM(name), viewer, ierr); CHKERRQ(ierr)
  call MatView(K, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'write_matrix'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine write_matrix
