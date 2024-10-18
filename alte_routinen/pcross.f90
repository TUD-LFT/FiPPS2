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
function pcross(vec1,vec2)
! =================================================================================================
!
!	Header:		function, return the cross product of two one dimensional arrays
!
!	Content:	function that computes the cross product of two one dimensional arrays
! 			with size (3,1)
!
!	Input:		vec1	- vector 1 in (vec1)x(vec2)
! 			vec2	- vector 2 in (vec1)x(vec2)
!
!	Output:		pcross	- cross product (vec1)x(vec2)
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
!	Revision
!
! =================================================================================================
!
! Use
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
! =================================================================================================
!
! Data types

! input

double precision, dimension(3), intent(in)	:: vec1, vec2

! output

double precision, dimension(3), intent(out)	:: pcross

! inner

double precision		:: size_vec
integer				:: size1, size2
integer				:: err_code=0
!
! =================================================================================================
!
! Calculation
!
! check for size (3,1) of arrays
!
if (size(vec1,1) == 3 .and. size(vec1,2) == 1 .or. size(vec2,1) == 3 .and. size(vec2,2) == 1) then

   ! compute cross product
   
   pcross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
   pcross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
   pcross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
   
else

   write(*,*) 'Error computing cross product, wrong array dimensions'
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
   write(*,*)  			   'pcross'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end function pcross
