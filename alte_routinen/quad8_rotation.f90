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
subroutine quad8_rotation (quad8_node_coords,mode,TRMat,quad8_node_coords_lo)
!
! =================================================================================================
!
!	Header:		subroutine for creating transformation matrix local -> global
!
!	Content:	Subroutine computes the transformation matrix from local to global coordinate
! 			system for quad-element based on nodal positions in global coordinate system
!
!	Input:		quad4_node_coords	- node-ids and nodal coordinates of quad4-element nodes
! 			ndof_			- number of degrees of freedom
!
!	Output:		TRLG	- transformation matrix from local to global coordinates
! 			Ael	- element surface
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
use netz_variablen
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
double precision, dimension(8,4), intent(in)	:: quad8_node_coords
character(2)					:: mode
!
! Output
!
double precision, dimension(48,48), intent(out)	:: TRmat
double precision, dimension(8,4),intent(out)	:: quad8_node_coords_lo
!
! Input + Output
!

!
! inner
!
double precision, dimension(3)			:: AC, BD
double precision, dimension(3)			:: e1, e2, e3, vec3
integer						:: imat, jmat
double precision, dimension(3,3)		:: BC

integer						:: kk,ll,mm,zz
integer						:: err_code=0
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
! Calculate vectors defining quadrilateral element
!
!AB(1:3) = quad8_node_coords(2,2:4)-quad8_node_coords(1,2:4)
AC(1:3) = quad8_node_coords(3,2:4)-quad8_node_coords(1,2:4)
!AD(1:3) = quad8_node_coords(4,2:4)-quad8_node_coords(1,2:4)

BD(1:3) = quad8_node_coords(4,2:4)-quad8_node_coords(2,2:4)

AC = AC/(sqrt(dot_product(AC,AC)))
BD = BD/(sqrt(dot_product(BD,BD)))

! Calculate element domain		 - not needed for calculation of stiffness matrix

!Ael = norm(pcross(AB,AD))

! Calculate unit-vectors

! nach Li.Z.X. A 9-node co-rotational quadrilateral shell element

e1 = AC-BD
! e1 = e1/norm(e1)
e1 = e1/(sqrt(dot_product(e1,e1)))

e2 = AC+BD
! e2 = e2/norm(e2)
e2 = e2/(sqrt(dot_product(e2,e2)))

!e3 = pcross(e1,e2)
e3(1) = (e1(2)*e2(3) - e1(3)*e2(2))
e3(2) = (e1(3)*e2(1) - e1(1)*e2(3))
e3(3) = (e1(1)*e2(2) - e1(2)*e2(1))

!
! Transformation submatrix
!
BC(1:3,1)=e1
BC(1:3,2)=e2
BC(1:3,3)=e3

!write(*,*) 'e1 :', e1
!write(*,*) 'e2 :', e2
!write(*,*) 'e3 :', e3
!write(*,*) 'BC :', BC
!stop

!
! Blow up to transformation matrix
!
do kk=1,(8*ndof/3)	! 8*3 Translations, 8*3 Rotations
   imat=(kk-1)*3
   jmat=(kk-1)*3
   do ll=1,3
      do mm=1,3
         TRMat(imat+ll,jmat+mm) = BC(ll,mm)
      end do
   end do
end do
!
! possible control: det(TRLG) != 1
!

! get nodal coordinates in local coordinate system
! defined via translation to coordinate system origin (A) -> AB, AC, AD
!  and rotation of direction vector global -> local

BC = transpose(BC)

! corner nodes

!quad8_node_coords_lo(1,1)   = 1.D0
quad8_node_coords_lo(1,1)   = quad8_node_coords(1,1)
quad8_node_coords_lo(1,2:4) = 0.D0

do kk=2,4
   
   ! get vectors from local coordinate system origin to corner nodes
   
   vec3 = quad8_node_coords(kk,2:4)-quad8_node_coords(1,2:4)
   
   ! rotate vector in local coordinate system
   
   vec3 = matmul(BC,vec3)
   
   ! write local coordinates
   
   quad8_node_coords_lo(kk,1)   = quad8_node_coords(kk,1)
   quad8_node_coords_lo(kk,2:4) = vec3
   
end do

! mid-side nodes

do kk=1,4
   
   zz=kk+1
   
   if (zz > 4) then
      zz = zz-4
   end if
      
   ! get vectors of quadrilateral sides in local coordinates
   
   vec3 = quad8_node_coords_lo(zz,2:4)-quad8_node_coords_lo(kk,2:4)
   
   ! write local coordinates
   
   quad8_node_coords_lo(kk+4,1)   = quad8_node_coords(kk+4,1)
   quad8_node_coords_lo(kk+4,2:4) = quad8_node_coords_lo(kk,2:4)+0.5*vec3
   
end do

! transformation matrix so far is for trafo local -> global

if (mode == 'lg') then

! do nothing

else if (mode == 'gl') then

   TRMat = transpose(TRMat)
   
else
   
   write(*,*) 'Error computing transformation matrix'
   write(*,*) 'Input .mode. falsely defined'
   err_code = 1
   goto 9999
   
end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)  			   'quad8_rotation'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if
return

end subroutine quad8_rotation
