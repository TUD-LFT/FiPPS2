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
subroutine lsolid20_extrapolate_gp(values_gp, values_cn)
! =================================================================================================
!
!   Header:     Exxtrapolate strains or stresses from 2x2 gauss points in a xi-eta-plane to the
!               corner nodes of the same plane
!
!   Content:    Input values are extrapolated from the 2x2 gauss points in a xi-eta-plane to the
!               corner nodes of the same plane; input values must be given with respect to the
!               element coordinate system;
!
!   Input:      values_gp: (6,4)-array containing the strains or stresses at the plane's 2x2 gauss
!                          points with respect to the element coordinate system
!
!   Output:     values_cn: (6,4)-array containing the strains or stresses at the plane's corner
!                          nodes with respect to the element coordinate system
!
!   Calls:      -
!
!   Called by:  lsolid20_strains_stresses
!
!   Author:     Florian Dexl
!               TU Dresden, Diplomarbeit
!          #LS20Marker
!
! =================================================================================================
!
! use
!
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
double precision, intent(in), dimension(6,4)    :: values_gp
!
! Output
!
double precision, intent(out), dimension(6,4)   :: values_cn
!
! Internal
!
integer                                         :: err_code=0
integer                                         :: ii
double precision, dimension(4,4)                :: EXT
double precision, parameter                     :: fac=1.d0/sqrt(3.d0)
!
! =================================================================================================
!
! Initialisation
!
 ! Matrix for extrapolation from integration points to element nodes
 EXT(1,1) = 1.d0+sqrt(3.d0)/2.d0; EXT(1,2) = -0.5d0
 EXT(1,3) = 1.d0-sqrt(3.d0)/2.d0; EXT(1,4) = -0.5d0
 EXT(2,1) = -0.5d0;               EXT(2,2) = 1.d0+sqrt(3.d0)/2.d0
 EXT(2,3) = -0.5d0;               EXT(2,4) = 1.d0-sqrt(3.d0)/2.d0
 EXT(3,1) = 1.d0-sqrt(3.d0)/2.d0; EXT(3,2) = -0.5d0
 EXT(3,3) = 1.d0+sqrt(3.d0)/2.d0; EXT(3,4) = -0.5d0
 EXT(4,1) = -0.5d0;               EXT(4,2) = 1.d0-sqrt(3.d0)/2.d0
 EXT(4,3) = -0.5d0;               EXT(4,4) = 1.d0+sqrt(3.d0)/2.d0
!
! =================================================================================================
!
! Calculation
!
 ! extrapolate stresses from gauss-points to element corner nodes
 do ii = 1,6
   values_cn(ii,:) = matmul(EXT,values_gp(ii,:))
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_extrapolate_gp'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_extrapolate_gp