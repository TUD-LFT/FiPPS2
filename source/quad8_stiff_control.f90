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
!> control subroutine for computing element stiffness matrix for quadratic 
!> eight-node serendipity shell/plate element and assemble it to the total
!> stiffness matrix
!
!> @details
!> subroutine calculates the element stiffness matrix for quadratic eight-node
!> quadrilateral shell/plate element and assembles it to the total stiffness matrix
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit, 08.12.2009
!               
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter,20.04.2017
!
!> $Id: quad8_stiff_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_stiff_control (fesim, KaaS)
! =================================================================================================
!
! Use
!
use konstanten
use fesimulation_typen
use globale_variablen
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
type(fe_simulation)                     :: fesim
!
! Output
!

!
! Input + Output
!
Mat                                     :: KaaS           !< AA-Part of of the global stiffness matrix (sparse - PETSC)
!
! inner
!
double precision, dimension(48,48)      :: quad8_K
double precision, dimension(8,3)        :: node_coords

double precision                        :: area, areaWeight
double precision                        :: E, nue, fac

integer                                 :: kk, jj, facI
integer                                 :: err_code=0

PetscMPIInt                             :: rank
PetscErrorCode                          :: ierr

double precision, dimension(6,6,1)      :: CTemp, TTemp
double precision, dimension(1)          :: lthTemp
integer, dimension(1)                   :: nopTemp
!
! =================================================================================================
!
! Initialisation
!
facI   = 1
!
! =================================================================================================
!
! Calculation
!
call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)


! Loop over all flat (reduced) quad8 elements
do jj=1,size(fesim%elemente%quad8s,1)

  if (mpi_quad8_proc_dist(jj) .NE. rank) cycle

   if (textoutput == .true. .and. jj == facI*100) then
      
      write(*,*) 'aktuelles Quad8-Element: ', facI*100
      
      facI = facI+1
      
   end if

   do kk = 1,8
     node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%quad8s(jj)%nids(kk))%coords(1:3)
   end do

   if (fesim%elemente%quad8s(jj)%propType == 1) then
   
    E   = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%intMat1ID)%ym
    nue = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%intMat1ID)%nu
      
    fac = E / (1.d0 - nue**2)
    
    CTemp = 0.d0    
    CTemp(1,1,1) = fac;          CTemp(1,2,1) = nue*fac;
    CTemp(2,1,1) = CTemp(1,2,1); CTemp(2,2,1) =     fac;
    CTemp(3,3,1) = 0.d0
    CTemp(4,4,1) = fac * (1.d0 - nue)/2.d0
    CTemp(5,5,1) = fac * (1.d0 - nue)/2.d0
    CTemp(6,6,1) = CTemp(5,5,1)
    
    TTemp = 0.d0
    do kk = 1, 6
      TTemp(kk,kk,1) = 1.d0
    end do
    
    lthTemp = fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%mt
    
    nopTemp = 2
   
    call quad8_stiff(1, &
                   & CTemp, &
                   & TTemp, &
                   & lthTemp, &
                   & fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%mt, &
                   & E, &
                   & node_coords, &
                   & nip_quad8, &
                   & nopTemp, &
                   & 93, &
                   & fesim%elemente%quad8s(jj)%area, &
                   & jj, &
                   & fesim%elemente%quad8s(jj)%eid, &
                   & quad8_K)
      
  else if (fesim%elemente%quad8s(jj)%propType == 2) then
   
    call quad8_stiff(fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%lay, &
                   & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%C(:,:,:), &
                   & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%T(:,:,:), &
                   & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%lth(:), &
                   & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%tth, &
                   & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%penaltyStiffness, &
                   & node_coords, &
                   & nip_quad8, &
                   & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%nop(:), &
                   & 91, &
                   & fesim%elemente%quad8s(jj)%area, &
                   & jj, &
                   & fesim%elemente%quad8s(jj)%eid, &
                   & quad8_K)

   else
      write(*,*) 'wrong property type for quad8-element',jj
      err_code=1
      goto 9999
   end if
   
   area = fesim%elemente%quad8s(jj)%area

   if (fesim%elemente%quad8s(jj)%propType == 1) then
      areaWeight = fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%areaWeight
      fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%weight = fesim%eigenschaften%pshells(fesim%elemente%quad8s(jj)%int_pid)%weight + area*areaWeight
   else if (fesim%elemente%quad8s(jj)%propType == 2) then
      areaWeight = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%areaWeight
      fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%weight = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(jj)%int_pid)%weight + area*areaWeight
   end if

   ! assemble element stiffness matrix in total stiffness matrix
   call assemble_aa (fesim, 8, fesim%elemente%quad8s(jj)%nids(1:8), quad8_K, KaaS)
   
!    call write_quad_matrix(quad8_K,48,"quad8_k             ","ANSYS")
     
end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_stiff_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_stiff_control
