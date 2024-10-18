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
!> $Id: sol1_solve.f90 381 2018-11-01 15:05:25Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 381 $
!> $Date: 2018-11-01 16:05:25 +0100 (Do, 01. Nov 2018) $
!
! =================================================================================================
subroutine sol1_calculate_tse(fesim, K, u, rank)
!
! Use
!
use konstanten
use fesimulation_typen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscksp.h"
  use petscksp
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  type(fe_simulation)                   :: fesim
  PetscMPIInt, intent(in)               :: rank
  Mat, intent(in)                       :: K
  Vec, intent(in)                       :: u
  Vec                                   :: h          ! Hilfsvektor als MPI-Sparse-Vektor
  PetscScalar                           :: w

  PetscErrorCode                        :: ierr

  integer                               :: err_code=0
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

 ! Berechnung der Gesamtdehnungsenergie (Total Strain Energy)
 if (rank .eq. 0) write(*,*) 'Berechne Gesamtdehnungsenergie (Total Strain Energy)'
 call VecDuplicate(u,h,ierr); CHKERRA(ierr)        ! Anlegen des Hilfsvektors  
 call MatMult(K,u,h,ierr); CHKERRA(ierr)
 call VecDot(u,h,w,ierr); CHKERRA(ierr)
 call VecDestroy(h,ierr); CHKERRA(ierr)            ! Zerstoeren des Hilfsvektors
 fesim%ergebnisse%tse = 0.5d0*w

 if (rank .eq. 0) write(*,*) 'TSE: ', fesim%ergebnisse%tse

!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol1_calculate_tse'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine sol1_calculate_tse
