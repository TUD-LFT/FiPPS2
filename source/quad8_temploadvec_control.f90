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
!> Adds temperature loadvector of all quad8 to the total loadvector
!
!> @details
!> Temperature load vector of each element is calculated by the corresponding
!> subroutine and then added to the total loadvector
!
!> @author
!> Andreas Hauffe, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 23.05.2017
!
!> $Id: quad8_temploadvec_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_temploadvec_control(fesim, lvec, pos, nodal_temperatures, quad8temperatures)
!
! use
!
use konstanten
use fesimulation_typen
!
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
 integer, intent(in)                             :: pos                             !< Position of loadvector in loadvector array (subcase)
 double precision, intent(in)                    :: nodal_temperatures(fesim%num_nodes)   !< nodal temperatures
 double precision, intent(in)                    :: quad8temperatures(size(fesim%elemente%quad8s,1))  !< elemental temperatures
!
! Input + Output
 double precision, intent(inout)                 :: lvec(1:fesim%num_dof,1:fesim%num_loadcases) !< Total loadvector for current load
!
! Internal
!
integer                                          :: ind_begin, lvec_begin
double precision                                 :: ncoords(8,3), ntemps(8), templvec(48)
double precision                                 :: transMat(3,3)

integer                                          :: err_code=0
integer                                          :: ii, jj, kk, nodeID
double precision, dimension(6,6,1)               :: CTemp, TTemp
double precision, dimension(1)                   :: lthTemp
integer, dimension(1)                            :: nopTemp
double precision, dimension(6,1)                 :: ath
double precision                                 :: E, nue, fac
!
! =================================================================================================
!
! Calculation
!
  do ii = 1,size(fesim%elemente%quad8s,1)
      
    ! Get node coordinates
    do jj = 1,8
      ncoords(jj,1:3) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(jj))%coords(1:3)
    end do
      
    if (fesim%is_quad8temp .eqv. .true.) then
      ! Get Temperature at elements
      ntemps(:)       = quad8temperatures(ii)
    else
      ! Get Temperature at element nodes
      ! Set zero, if no temperature is prescribed
      ! Get coordinates of element nodes
      ntemps = 0.d0
      do jj = 1,8
        ntemps(jj)    = nodal_temperatures(fesim%elemente%quad8s(ii)%nids(jj))
      end do
    end if
    
    if (fesim%elemente%quad8s(ii)%propType == 1) then
    
      E   = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(ii)%int_pid)%intMat1ID)%ym
      nue = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(ii)%int_pid)%intMat1ID)%nu
        
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
      
      lthTemp = fesim%eigenschaften%pshells(fesim%elemente%quad8s(ii)%int_pid)%mt
      
      nopTemp = 2
      
      ath(1:3,1) = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(ii)%int_pid)%intMat1ID)%ath
      ath(4:6,1) = 0.d0
    
      call quad8s_temploadvec(1, &
                    & CTemp, &
                    & TTemp, &
                    & ath, &
                    & lthTemp, &
                    & fesim%eigenschaften%pshells(fesim%elemente%quad8s(ii)%int_pid)%mt, &
                    & ncoords, &
                    & nip_quad8, &
                    & nopTemp, &
                    & 93, &
                    & ntemps, &
                    & templvec)
        
    else if (fesim%elemente%quad8s(ii)%propType == 2) then

      call quad8s_temploadvec(fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%lay, &
                    & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%C(:,:,:), &
                    & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%T(:,:,:), &
                    & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%ath(:,:), &
                    & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%lth(:), &
                    & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%tth, &
                    & ncoords, &
                    & nip_quad8, &
                    & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(ii)%int_pid)%nop(:), &
                    & 91, &
                    & ntemps, &
                    & templvec)

    else
        write(*,*) 'wrong property type for quad8-element',ii
        err_code=1
        goto 9999
    end if
    
    do jj = 1,8
      nodeID     = fesim%elemente%quad8s(ii)%nids(jj)
      
      ind_begin  = (jj-1)*6+1
      
      lvec_begin = (nodeID-1)*ndof+1

      if (fesim%knoten%nodes(nodeID)%cid .NE. 0) then
        transMat(1:3,1:3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nodeID)%cid)%transMat(1:3,1:3)
        transMat = transpose(transMat)
        templvec(ind_begin  :ind_begin+2) = matmul(transMat,templvec(ind_begin  :ind_begin+2))
        templvec(ind_begin+3:ind_begin+5) = matmul(transMat,templvec(ind_begin+3:ind_begin+5))
      end if

      lvec(lvec_begin  :lvec_begin+2,pos) = lvec(lvec_begin  :lvec_begin+2,pos)+templvec(ind_begin  :ind_begin+2)
      lvec(lvec_begin+3:lvec_begin+5,pos) = lvec(lvec_begin+3:lvec_begin+5,pos)+templvec(ind_begin+3:ind_begin+5)
    end do
    
  end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_temploadvec_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_temploadvec_control
