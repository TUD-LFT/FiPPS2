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
!
!> @details
!
!> @author 
!
!> $Id: init_vectors.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine init_vectors(dim_dof, FaS, UaS)

#include "petsc/finclude/petscvec.h"
  use petscvec
!

!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  integer, intent(in)           :: dim_dof
  
  PetscErrorCode                :: ierr

  Vec                           :: UaS          ! Verschiebungsvektor als MPI-Sparse-Vektor
  Vec                           :: FaS          ! Kraftvektor ohne gesperrte Freiheiten als Sparse-Vektor
  
  integer                       :: err_code=0
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
   ! Anlegen von Fa
   call VecCreate(PETSC_COMM_WORLD,FaS,ierr); CHKERRQ(ierr)
   call VecSetSizes(FaS,PETSC_DECIDE,dim_dof,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(FaS,ierr); CHKERRQ(ierr)
 
   ! Anlegen von Ua
   call VecCreate(PETSC_COMM_WORLD,UaS,ierr); CHKERRQ(ierr)
   call VecSetSizes(UaS,PETSC_DECIDE,dim_dof,ierr); CHKERRQ(ierr)
   call VecSetFromOptions(UaS,ierr); CHKERRQ(ierr)
 
!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'get_Utot'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if
  
end subroutine init_vectors
