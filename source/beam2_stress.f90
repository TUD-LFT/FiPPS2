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
!> subroutine for stress calculation of 2-node beam element
!
!> @details
!> subroutine calculates the stresses of a 2-node straight beam element
!> with constant cross-section properties
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 20.09.2010
!
!> $Id: beam2_stress.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_stress(fesim, dispg, nodal_temperatures, beam2temperatures, stress, thermal_force, thermal_stress)
! =================================================================================================
!
!       Header:         
!
!       Content:        
!
!       Input:          dispg                   - displacements in global coordinates
!
!       Output:         stress                  - stresses for each beam element
!
!       Internal:       
!
!       Calls:          beam2_get_properties
!
!       Called by:      sol1_output
!
!       Author:         Andreas Hauffe                  20.09.2010
!                       TU Dresden, WiMi
!
!       Revision:       
!
! =================================================================================================
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
! =================================================================================================
!
! Data types
  type(fe_simulation)                                      :: fesim                  !< FE-simulation data
!
! Input
!  
double precision, dimension(:), intent(in)                 :: dispg(fesim%num_dof)   !< displacement in global coordinates
double precision, dimension(fesim%num_nodes), intent(in)   :: nodal_temperatures     !< nodal temperatures
double precision, dimension(size(fesim%elemente%beam2s,1)), intent(in):: beam2temperatures !< elemental temperatures
!
! Output
!
double precision, dimension(:), intent(out)                :: stress(fesim%num_elements,7) !< stresses for each beam element
double precision, dimension(:), intent(out)                :: thermal_force(fesim%num_elements)
double precision, dimension(:), intent(out)                :: thermal_stress(fesim%num_elements)
! Input + Output
!

!
! Inner
!
integer                             :: err_code=0
double precision, dimension(2,3)    :: beam2_node_coords
integer, dimension(2)               :: beam2_node_ids
double precision, dimension(12)     :: disp_beam2
double precision                    :: AA,I11,I22,I12,It,Emat,Gmat,alpha
double precision, dimension(3)      :: vv
double precision, dimension(12,12)  :: TRMat
double precision                    :: le
double precision                    :: My1,My2,Mz1,Mz2,Fx,Fy,Fz
double precision                    :: str_1, str_2, ty, tz
double precision                    :: btemp, temp_force

integer                             :: pos,ii,jj,kk

! Loop over all 2-node 3D beam elements
do ii=1, size(fesim%elemente%beam2s,1)
   
   ! get node locations
   
   do jj=1,2
      
      beam2_node_ids(jj) = fesim%elemente%beam2s(ii)%nids(jj)
      beam2_node_coords(jj,1:3) = fesim%knoten%nodes(fesim%elemente%beam2s(ii)%nids(jj))%coords(1:3)
      
   end do
   
   ! get global displacement vector for element nodes

   do jj=1,2
      
      pos = (fesim%elemente%beam2s(ii)%nids(jj)-1)*ndof+1
      
      do kk=1,ndof
        disp_beam2(kk+(jj-1)*ndof) = dispg(pos+kk-1)
      end do
      
   end do
      
   ! get beam2-element material stiffness and section properties
   
   call beam2_get_properties(fesim,ii,AA,I11,I22,I12,It,Emat,Gmat)
   alpha = fesim%materialien%mat1s(fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%intMat1ID)%ath
   
   do jj=1,3
     vv(jj)=fesim%elemente%beam2s(ii)%xi(jj)
   end do

   ! get transformation matrix from global to local coordinate system
   
   call beam2_rotation (beam2_node_coords,vv,'gl',TRMat,le)
   
   ! transform displacement vector to local coordinates
   
   disp_beam2=matmul(TRMat,disp_beam2)

   ! calculate element moments and forces
   
!***********************************
!
!  Die Momente müssen noch um den Winkel theta rotiert werden.
!
!***********************************
   My1 = (Emat*I22/(le**2))*(6.D0*(disp_beam2(9)-disp_beam2(3))+4.D0*le*disp_beam2(5)+2.D0*le*disp_beam2(11))
   My2 = (Emat*I22/(le**2))*(6.D0*(disp_beam2(9)-disp_beam2(3))+2.D0*le*disp_beam2(5)+4.D0*le*disp_beam2(11))
   
   Mz1 = (Emat*I11/(le**2))*(-6.D0*(disp_beam2(8)-disp_beam2(2))+4.D0*le*disp_beam2(6)+2.D0*le*disp_beam2(12))
   Mz2 = (Emat*I11/(le**2))*(-6.D0*(disp_beam2(8)-disp_beam2(2))+2.D0*le*disp_beam2(6)+4.D0*le*disp_beam2(12))
   
    ! Get Temperature at element nodes
    ! Set zero, if no temperature is prescribed
    ! Get coordinates of element nodes
    btemp = 0.d0
    if (fesim%is_beam2temp == .true.) then
      btemp = beam2temperatures(ii)
    else if (fesim%is_temperature == .true.) then
      do jj = 1,2
        btemp = btemp + nodal_temperatures(fesim%elemente%beam2s(ii)%nids(jj))
      end do
      btemp = btemp/2.d0
    end if

   temp_force = Emat*AA*btemp*alpha

   thermal_force(fesim%elemente%beam2s(ii)%eid) = temp_force
   thermal_stress(fesim%elemente%beam2s(ii)%eid) = temp_force/AA
   
   Fx = (Emat*AA/le)*(disp_beam2(7)-disp_beam2(1)) - temp_force
   Fy = (-1.D0/le)*(Mz2+Mz1)
   Fz = ( 1.D0/le)*(My2+My1)
   
   stress(fesim%elemente%beam2s(ii)%eid,1) = Fx/AA
   
   tz = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%t1
   ty = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%t2
   
   str_1 = ABS(My1 * tz / (2.d0 * I22))
   str_2 = ABS(My2 * tz / (2.d0 * I22))
   stress(fesim%elemente%beam2s(ii)%eid,2) = MAX(str_1, str_2)
   
   str_1 = ABS(Mz1 * ty / (2.d0 * I11))
   str_2 = ABS(Mz2 * ty / (2.d0 * I11))
   stress(fesim%elemente%beam2s(ii)%eid,3) = MAX(str_1, str_2)
   
   stress(fesim%elemente%beam2s(ii)%eid,4) = stress(fesim%elemente%beam2s(ii)%eid,1) + stress(fesim%elemente%beam2s(ii)%eid,2) + stress(fesim%elemente%beam2s(ii)%eid,3)
   
   stress(fesim%elemente%beam2s(ii)%eid,5) = stress(fesim%elemente%beam2s(ii)%eid,1) - stress(fesim%elemente%beam2s(ii)%eid,2) - stress(fesim%elemente%beam2s(ii)%eid,3)
   
   stress(fesim%elemente%beam2s(ii)%eid,6) = Fx
   
   stress(fesim%elemente%beam2s(ii)%eid,7) = MAX(MAX(ABS(My1), ABS(My2)), MAX(ABS(Mz1), ABS(Mz2)))
   
end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_stress'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine beam2_stress
