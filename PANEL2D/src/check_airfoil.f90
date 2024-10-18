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
SUBROUTINE check_airfoil()

use functions
use konst_var
use netz_variablen

implicit none

! Declaration of variables
integer                                           :: ii
logical                                           :: closed_te
logical                                           :: double_points

double precision, parameter                       :: eps=1.d-10

! Calculation

! Check for closed trailing edge
closed_te = (abs_vector2(nodes(1,:) - nodes(n_nodes,:)) .lt. eps)

if (closed_te .neqv. .true.) then
    write(*,*) '  ERROR in airfoil input data!'
    write(*,*) '  Trailing edge must be closed.'
    STOP
end if

! Check for double defined points
double_points = .false.
do ii = 2, n_nodes
    if (abs_vector2(nodes(ii,:) - nodes(ii-1,:)) .lt. eps) then
        double_points = .true.
        exit
    end if
end do

if (double_points .neqv. .false.) then
    write(*,*) '  ERROR in airfoil input data!'
    write(*,*) '  Points double defined.'
    STOP
end if

END SUBROUTINE check_airfoil
