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
!
!> @details
!
!> @author 
!
!> $Id: vec_func.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module vec_func_own


contains

function pcross(vec1,vec2)
!
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
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!

! input

double precision, dimension(3), intent(in)	:: vec1, vec2

! output

double precision, dimension(3)           	:: pcross

!
! =================================================================================================
!
! Calculation
!
! check for size (3,1) of arrays
!

   ! compute cross product
   
   pcross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
   pcross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
   pcross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

end function pcross

end module vec_func_own
