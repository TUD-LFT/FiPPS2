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
!> Control subroutine for computing elemental strain energyies for quadratic twenty-node
!> layered solid element lsolid20
!
!> @details
!> Subroutine controls the ca
!> Elementwise calculation of the total strain energy (tse) for lsolid20 elements.
!
!> @author   Florian Dexl, TU Dresden, wiss. Mitarbeiter, 16.10.2023
!
!> $Id: quad8_geostiff_control.f90 421 2020-03-02 15:22:57Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 421 $
!> $Date: 2020-03-02 16:22:57 +0100 (Mo, 02. Mär 2020) $
!
!
! =================================================================================================
subroutine lsolid20_strainEnergy_control (fesim, dispg, scloop, tse)
! =================================================================================================
! Use
!
use konstanten
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
type(fe_simulation)                                                 :: fesim
double precision, dimension(:), intent(in)                          :: dispg(fesim%num_dof)         !< (num_dof)-Array containing the displacements in the nodal coordinate system
integer, intent(in)                                                 :: scloop
!
! Output
!
double precision, dimension(fesim%num_elements), intent(out)        :: tse            !< total strain energy of elements
!
! inner
!
double precision, dimension(20*3)                                   :: displ
double precision                                                    :: lsolid20_tse

integer                                                             :: ind, kk, ll, jj
integer                                                             :: nip, nlay, cid
integer                                                             :: node_id
integer                                                             :: err_code=0
double precision, dimension(20,3)                                   :: node_coords
double precision                                                    :: tth
double precision, dimension(:,:,:),allocatable                      :: Clay
double precision, dimension(:),allocatable                          :: lth
integer, dimension(:), allocatable                                  :: nop
!
! =================================================================================================
!
! Initialisation
!
 nip = 2    ! number of integration points in plane
!
! =================================================================================================
!
! Calculation
!

! Loop over all layered solid elements
 do jj = 1,size(fesim%elemente%lsolid20s,1)

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
     node_id        = fesim%elemente%lsolid20s(jj)%nids(kk)
     node_coords(kk,1:3) = fesim%knoten%nodes(node_id)%coords(1:3)
     ind = (node_id-1)*ndof
     do ll=1,3
       displ(ll+(kk-1)*3) = dispg(ind+ll)
     end do
   end do
   
   ! calculate element strain energy 
   call lsolid20_strainEnergy(fesim,nlay,Clay,lth,tth,node_coords,nip,nop,cid,displ,jj,scloop,lsolid20_tse)
   
   deallocate(Clay,lth,nop)

   tse(fesim%elemente%lsolid20s(jj)%eid) = lsolid20_tse

end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_strainEnergy_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_strainEnergy_control
