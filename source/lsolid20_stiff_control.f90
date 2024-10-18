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
!> layered 20-node solid element and assemble it to the total
!> stiffness matrix
!
!> @details
!> subroutine calculates the element stiffness matrix for quadratic layered 20-node
!> solid element and assembles it to the total stiffness matrix
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_stiff_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_stiff_control (fesim, KaaS)

! =================================================================================================
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
type(fe_simulation), intent(in)         :: fesim
!
! Output
!

!
! Input + Output
!
Mat                 :: KaaS !< KaaS
!
! Internal
!
integer, dimension(20)              :: node_ids
double precision, dimension(20,3)   :: node_coords
double precision, dimension(60,60)  :: lsolid20_K
double precision, dimension(120,120):: K
double precision, allocatable       :: Clay(:,:,:), lth(:)
integer, allocatable, dimension(:)  :: nop

double precision                    :: tth

integer                             :: jj, kk, facI, nlay, nip, cid
integer                             :: err_code=0

PetscMPIInt                         :: rank
PetscErrorCode                      :: ierr
!character(20)                       :: name
!
! =================================================================================================
!
! Initialisation
!
 facI   = 1
 nip = 2    ! number of integration points in plane
!
! =================================================================================================
!
! Calculation
!
 call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

! Loop over all layered solid elements
 do jj = 1,size(fesim%elemente%lsolid20s,1)
   if (mpi_lsolid20_proc_dist(jj) .NE. rank) cycle

   if (textoutput == .true. .and. jj == facI*100) then
    
     write(*,*) 'aktuelles Lsolid20-Element: ', facI*100
    
     facI = facI+1
    
   end if

   ! get number of layers for current element
   nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%lay
   
   ! allocate arrays dependent on layer number
   allocate(Clay(6,6,nlay))
   allocate(lth(nlay))
   allocate(nop(nlay))
   
   ! get stiffness matrices for layers
   Clay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%C(:,:,:)
   ! get thickness for layers
   lth  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%lth(:)
   ! get total thickness of laminat
   tth  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%tth
   ! get number of integration points over thickness for each layer
   nop  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%nop(:)
   ! get ID of coordinate system for orientation of element coosy
   cid  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%cid

   do kk = 1,20
     node_ids(kk)        = fesim%elemente%lsolid20s(jj)%nids(kk)
     node_coords(kk,1:3) = fesim%knoten%nodes(node_ids(kk))%coords(1:3)
   end do
   
   ! calculate element stiffness matrix  
   call lsolid20_stiff(fesim,nlay,Clay,lth,tth,node_coords,nip,nop,cid,jj,fesim%elemente%lsolid20s(jj)%eid,lsolid20_K)

   ! Blow up stiffness matrix to dof = 6
   ! Zeros are added at entries of rotational dofs
   call blow_up_stiff(20, lsolid20_K, K)
   
   call assemble_aa (fesim,20,node_ids(:),K,KaaS)
   
   deallocate(Clay,lth,nop)

!    ! write matrix to standard out
    !write(name,'(A3,I17)') 'Kel',jj
    !call write_quad_matrix(lsolid20_K, 60, name, 'FiPPS')
    !call write_quad_matrix(lsolid20_K, 60, 'Kel                 ', 'FiPPS')
   
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_stiff_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_stiff_control
