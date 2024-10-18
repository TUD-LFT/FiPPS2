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
!> Adds temperature loadvector of all lSolid20 to the total loadvector
!
!> @details
!> Temperature load vector of each element is calculated by the corresponding
!> subroutine and then added to the total loadvector
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_temploadvec_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_temploadvec_control(fesim, lvec, pos, nodal_temperatures, lsolid20temperatures)
! =================================================================================================
! use
!
use konstanten
use fesimulation_typen
! =================================================================================================
!
  implicit none
! =================================================================================================
!
! Data types
!
! Input
!
 type(fe_simulation)                             :: fesim
 integer, intent(in)                             :: pos !< Position of loadvector in loadvector array
 double precision, intent(in)                    :: nodal_temperatures(fesim%num_nodes) !< nodal temperatures
 double precision, intent(in)                    :: lsolid20temperatures(size(fesim%elemente%lsolid20s,1))  !< elemental temperatures
!
! Input + Output
 double precision, intent(inout)                 :: lvec(1:fesim%num_dof,1:fesim%num_loadcases) !< Total loadvector for current load
!
! Internal
!
integer                                          :: nlay, nip, cid
integer                                          :: ind_begin, ind_end, lvec_begin, lvec_end
integer, allocatable, dimension(:)               :: nop
double precision, allocatable, dimension(:)      :: lth
double precision, allocatable, dimension(:,:)    :: ath_lay
double precision, allocatable, dimension(:,:,:)  :: Clay
double precision                                 :: tth
double precision                                 :: ncoords(20,3), ntemps(20), templvec(60)
double precision                                 :: transMat(3,3)

integer                                          :: err_code=0
integer                                          :: ii, jj, nodeID
!
! =================================================================================================
!
! Initialisation
!
  nip = 2    ! number of integration points in plane
! =================================================================================================
!
! Calculation
!
  do ii = 1,size(fesim%elemente%lsolid20s,1)
  
    ! Get node coordinates
    do jj = 1,20
      ncoords(jj,1:3) = fesim%knoten%nodes(fesim%elemente%lsolid20s(ii)%nids(jj))%coords(1:3)
    end do
  
    if (fesim%is_lsolid20temp .eqv. .true.) then
      ! Get Temperature at elements
      ntemps(:)       = lsolid20temperatures(ii)
    else if (fesim%is_temperature .eqv. .true.) then
      ! Get Temperature at element nodes
      ! Set zero, if no temperature is prescribed
      ! Get coordinates of element nodes
      ntemps = 0.d0
      do jj = 1,20
        ntemps(jj)    = nodal_temperatures(fesim%elemente%lsolid20s(ii)%nids(jj))
      end do
    end if
  
    ! get number of layers for current element
    nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%lay
    
    ! allocate arrays dependent on layer number
    allocate(Clay(6,6,nlay))
    allocate(ath_lay(6,nlay))
    allocate(lth(nlay))
    allocate(nop(nlay))
    
    ! get stiffness matrices for layers
    Clay    = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%C(:,:,:)
    ! get thermal expansion coefficients for layers
    ath_lay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%ath(:,:)
    ! get thickness for layers
    lth     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%lth(:)
    ! get total thickness of laminat
    tth     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%tth
    ! get number of integration points over thickness for each layer
    nop     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%nop(:)
    ! get ID of coordinate system for orientation of element coosy
    cid     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%cid

    call lsolid20_temploadvec(fesim,nlay,Clay,ath_lay,lth,tth,ncoords,nip,nop,cid,ntemps,templvec)
    
    do jj = 1,20
      nodeID     = fesim%elemente%lsolid20s(ii)%nids(jj)
      
      ind_begin  = (jj-1)*3+1
      ind_end    = ind_begin+2
      
      lvec_begin = (nodeID-1)*ndof+1
      lvec_end   = lvec_begin+2

      if (fesim%knoten%nodes(nodeID)%cid .NE. 0) then
        transMat(1:3,1:3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nodeID)%cid)%transMat(1:3,1:3)
        transMat = transpose(transMat)
        templvec(ind_begin:ind_end) = matmul(transMat,templvec(ind_begin:ind_end))
      end if

      lvec(lvec_begin:lvec_end,pos) = lvec(lvec_begin:lvec_end,pos)+templvec(ind_begin:ind_end)
    end do
    
    deallocate(Clay,ath_lay,lth,nop)
    
  end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_temploadvec_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_temploadvec_control
