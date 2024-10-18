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
subroutine testroutine (bnode,mode,BC)

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
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
! =================================================================================================
!

use globale_variablen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
double precision, dimension(3,4), intent(in)	:: bnode
character(2)					:: mode
!
! Output
!
double precision, dimension(3,3), intent(out)	:: BC
!
! Inner
!
double precision, dimension(3)			:: AB,AC,e1,e2,e3,pcABAC
!
! Input
!
!
! =================================================================================================
!
! Initialisation
!
!
! =================================================================================================
!
! Calculation
!

AB = bnode(2,2:4)-bnode(1,2:4)
AC = bnode(3,2:4)-bnode(1,2:4)

!
! Calc auxiliary values
!
pcABAC(1) = AB(2)*AC(3) - AB(3)*AC(2)
pcABAC(2) = AB(3)*AC(1) - AB(1)*AC(3)
pcABAC(3) = AB(1)*AC(2) - AB(2)*AC(1)

! Calc unit vectors

e3 = AB/(sqrt(dot_product(AB,AB)))
e2 = pcABAC/sqrt(dot_product(pcABAC,pcABAC))

e1(1) = -(e2(2)*e3(3) - e2(3)*e3(2))
e1(2) = -(e2(3)*e3(1) - e2(1)*e3(3))
e1(3) = -(e2(1)*e3(2) - e2(2)*e3(1))

! Build transformation matrix local -> global

BC(1:3,1) = e1
BC(1:3,2) = e2
BC(1:3,3) = e3

if (mode == 'lg') then

! do nothing

else if (mode == 'gl') then

   BC = transpose(BC)
   
else
   
   write(*,*) 'Error computing transformation matrix'
   write(*,*) 'Input .mode. falsely defined'
   goto 9999
   
end if

!
! =================================================================================================
!
! Error handling
!
9999 continue
write(*,*) 'An error occured, exit program '
return

end subroutine testroutine
