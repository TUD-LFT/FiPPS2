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
subroutine rotation_coord(bnode,mode,BC)

! =================================================================================================
!
!	Header:		subroutine for computing the transformation matrix from  a cartesian to the
!			global coordinate system
!
!	Content:	subroutine calculates the transformations matrix from an arbitrary cartesian
!			coordinate system, defined by 3 nodes, to the global cartesian coordinate
!			system
!			point 1:	origin
!			point 2:	point on z-axis of coordinate system
!			point 3:	point in x-z-plane of coordiante system
!
!	Input:		
!
!	Output:		
!
!	Internal:	
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
!	Revision:	
!
! =================================================================================================
!
! Use
!
!use globale_variablen
!use vec_func_own
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
! Input
!
double precision, dimension(3,4), intent(in)	:: bnode
character(2), intent(in)			:: mode
!
! Output
!
double precision, dimension(3,3), intent(out)	:: BC
!
! Input + Output
!
!
! inner
!
double precision, dimension(3)			:: AB,AC,e1,e2,e3,pcABAC
!!double precision, dimension(3)			:: pcross
!!
!! =================================================================================================
!!
!! Initialisation
!!
!!
!! =================================================================================================
!!
!! Calculation
!!
!AB(1:3) = bnode(2,2:4)-bnode(1,2:4)
!AC = bnode(3,2:4)-bnode(1,2:4)
!
!! Calc auxiliary values
!
!!pcABAC = pcross(AB,AC)
!
!! Calc unit vectors
!
!e3 = AB/(sqrt(dot_product(AB,AB)))
!e2 = pcABAC/sqrt(dot_product(pcABAC,pcABAC))
!!e1 = -pcross(e2,e3)
!
!! Build transformation matrix
!
BC=1.d0
!BC(1:3,1) = e1
!BC(1:3,2) = e2
!BC(1:3,3) = e3
!
!! transformation matrix so far is for trafo lokal -> global
!
!if (mode == 'lg') then
!
!! do nothing
!
!else if (mode == 'gl') then
!
!   BC = transpose(BC)
!   
!else
!   
!   write(*,*) 'Error computing transformation matrix'
!   write(*,*) 'Input .mode. falsely defined'
!   goto 9999
!   
!end if
!
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)  			   'rotation_coord'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine rotation_coord
