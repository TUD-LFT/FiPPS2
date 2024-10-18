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
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter
!
!> $Id: beam2_temploadvec_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_temploadvec_control(fesim, lvec, pos, nodal_temperatures, beam2temperatures)
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
 type(fe_simulation)                             :: fesim                                    !< FE-simulation data
 integer, intent(in)                             :: pos                                      !< position of loadvector in loadvector array
 double precision, intent(in)                    :: nodal_temperatures(fesim%num_nodes)      !< nodal temperatures
 double precision, intent(in)                    :: beam2temperatures(size(fesim%elemente%beam2s,1))  !< elemental temperatures
!
! Input + Output
 double precision, intent(inout)                 :: lvec(1:fesim%num_dof,1:fesim%num_loadcases) !< total loadvector for current load 
!
! Internal
!

double precision                                 :: AA,I11,I22,I12,It,Emat,Gmat,alpha
integer                                          :: ind_begin, lvec_begin
double precision                                 :: etemp, templvec(12)
double precision                                 :: transMat(3,3)
double precision, dimension(12,12)               :: TRMat
double precision, dimension(2,3)                 :: beam2_node_coords
double precision, dimension(3)                   :: vv
double precision                                 :: le

integer                                          :: err_code=0
integer                                          :: ii, jj, nodeID
!
! =================================================================================================
!
! Initialisation
!
! =================================================================================================
!
! Calculation
!
  do ii = 1,size(fesim%elemente%beam2s,1)
  
    call beam2_get_properties(fesim,ii,AA,I11,I22,I12,It,Emat,Gmat)
    
    alpha = fesim%materialien%mat1s(fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%intMat1ID)%ath

    if (fesim%is_beam2temp .eqv. .true.) then
      ! Get Temperature at elements
      etemp = beam2temperatures(ii)
    else    
      ! Get Temperature at element nodes
      ! Set zero, if no temperature is prescribed
      ! Get coordinates of element nodes
      etemp = 0.d0
      do jj = 1,2
        etemp = etemp + nodal_temperatures(fesim%elemente%beam2s(ii)%nids(jj))
      end do
      etemp = etemp/2.d0
    end if

    vv(1:3)=fesim%elemente%beam2s(ii)%xi(1:3)
    
    do jj=1,2
      beam2_node_coords(jj,1:3) = fesim%knoten%nodes(fesim%elemente%beam2s(ii)%nids(jj))%coords(1:3)
    end do

   ! get transformation matrix from local to global coordinate system
  
   ! call rotation_beam2 (beam2_node_coords,(fesim%elemente%beam2s(ii)%xi(jj), jj=1,3),'lg',TRMat,le)
    call beam2_rotation (beam2_node_coords,vv,'lg',TRMat,le)
    
    templvec = 0.d0
    
    templvec(1) = -Emat*AA*etemp*alpha
    templvec(7) = Emat*AA*etemp*alpha
    
    templvec = matmul(TRMat,templvec)
    
          
    do jj = 1,2
      nodeID     = fesim%elemente%beam2s(ii)%nids(jj)
      
      ind_begin  = (jj-1)*6+1
      
      lvec_begin = (nodeID-1)*ndof+1

      if (fesim%knoten%nodes(nodeID)%cid .NE. 0) then
        transMat = transpose(fesim%koordinatensysteme%coords(fesim%knoten%nodes(nodeID)%cid)%transMat)
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
   write(*,*)                      'beam2_temploadvec_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine beam2_temploadvec_control
