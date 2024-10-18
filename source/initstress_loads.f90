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
!> Adds loadvector of initial stresses to the total loadvector
!
!> @details
!> The equivalent load vector of initial stresses is calculated by the corresponding
!> subroutine and then added to the total loadvector
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: initstress_loads.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine initstress_loads(fesim, lvec, pos)
! =================================================================================================
! use
!
use konstanten
use fesimulation_typen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Include
!
!
! =================================================================================================
!
! Data types
!
! Input
!
 type(fe_simulation)                             :: fesim
 integer, intent(in)                             :: pos !< Position of loadvector in loadvector array
!
! Input + Output
 double precision, intent(inout)                 :: lvec(1:fesim%num_dof,1:fesim%num_loadcases) !< Total loadvector for current load
!
! Internal
!
integer                                          :: nlay, nip
integer                                          :: ind_begin, ind_end, lvec_begin, lvec_end
integer, allocatable, dimension(:)               :: nop
double precision, allocatable, dimension(:)      :: lth
double precision                                 :: tth
double precision                                 :: ncoords(20,3), initstrvec(60)
double precision                                 :: transMat(3,3)

integer                                          :: err_code=0
integer                                          :: ii, jj, nodeID
!
! =================================================================================================
!
! Initialisation
 nip = 2    ! number of integration points in plane
!
! =================================================================================================
!
! Calculation
!
 ! Loop over Lsolid20 elements
 do ii = 1,size(fesim%elemente%lsolid20s,1)
 
   ! Get coordinates of element nodes
   do jj = 1,20
     ncoords(jj,1:3) = fesim%knoten%nodes(fesim%elemente%lsolid20s(ii)%nids(jj))%coords(1:3)
   end do
 
   ! get number of layers for current element
   nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%lay
   
   ! allocate arrays dependent on layer number
   allocate(lth(nlay))
   allocate(nop(nlay))

   ! get thickness for layers
   lth     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%lth(:)
   ! get total thickness of laminat
   tth     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%tth
   ! get number of integration points over thickness for each layer
   nop     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%nop(:)

   call lsolid20_initstressvec(fesim,ii,nlay,lth,tth,ncoords,nip,nop,initstrvec)
  
   do jj = 1,20
     nodeID     = fesim%elemente%lsolid20s(ii)%nids(jj)
     
     ind_begin  = (jj-1)*3+1
     ind_end    = ind_begin+2
     
     lvec_begin = (nodeID-1)*ndof+1
     lvec_end   = lvec_begin+2

     if (fesim%knoten%nodes(nodeID)%cid .NE. 0) then
       transMat(1:3,1:3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nodeID)%cid)%transMat(1:3,1:3)
       transMat = transpose(transMat)
       initstrvec(ind_begin:ind_end) = matmul(transMat,initstrvec(ind_begin:ind_end))
     end if
     
     lvec(lvec_begin:lvec_end,pos) = lvec(lvec_begin:lvec_end,pos)-initstrvec(ind_begin:ind_end)
   end do
   
   deallocate(lth,nop)
  
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'initstress_loads'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine initstress_loads
