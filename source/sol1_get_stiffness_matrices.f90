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
!> control subroutine for obtaining the global elastic stiffness matrix from
!> the different elements and their element stiffness matrices
!
!> @details
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit, 03.12.2009
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 06.05.2010
!
!> $Id: sol1_get_stiffness_matrices.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_get_stiffness_matrices(fesim,KaaS)
!
  use fesimulation_typen
!
! =================================================================================================
!
! Include
!
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
! Input
!  
  type(fe_simulation), intent(in)               :: fesim
!
! Inner
!
  integer                                       :: err_code=0
!
! Output
! 
  Mat                                           :: KaaS         ! Steifigkeitsmatrix ohne gesperrte Freiheiten als Sparse-Matrix

  PetscMPIInt                                   :: rank
  PetscErrorCode                                :: ierr
!
! =================================================================================================
!
! Initialisation
!
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
!
! =================================================================================================
!
! Calculation
!
  if (fesim%is_quad8 == .true.) then
    call quad8_stiff_control (fesim, KaaS)
  end if
  
  if (fesim%is_lsolid20 == .true.) then
    call lsolid20_stiff_control(fesim, KaaS)
  end if

  if (rank == 0) then  
    if (fesim%is_beam2 == .true.) then
      call beam2_stiff_control (fesim, KaaS)
    end if
  end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_get_stiffness_matrices'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_get_stiffness_matrices
