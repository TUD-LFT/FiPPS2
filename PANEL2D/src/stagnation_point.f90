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
SUBROUTINE stagnation_point(vp, ispu, ispl)

! Load modules
use konst_var
use netz_variablen
use functions

implicit none

! Declaration of variables
intent(in)                                        :: vp
intent(out)                                       :: ispu, ispl
double precision                                  :: xmin
double precision, dimension(n_nodes)              :: vp

logical                                           :: stagnation_found

integer                                           :: ispu, ispl
integer                                           :: ii

! Calculation

! Find stagnation point next to leading edge
stagnation_found=.false.
xmin = HUGE(1.d0)
do ii = 1, n_nodes-1
    if ((vp(ii) .eq. 0.d0) .and. (nodes(ii,1) .lt. xmin)) then
        ispu = ii
        ispl = ii
        xmin = nodes(ii,1)
        stagnation_found = .true.
    else if ((vp(ii) .lt. 0.d0) .and. (vp(ii+1) .gt. 0.d0) .and. (nodes(ii,1) .lt. xmin)) then
        ispu = ii+1
        ispl = ii
        xmin = nodes(ii,1)
        stagnation_found = .true.
    end if
end do

if (stagnation_found .EQV. .false.) then
  write(*,*) 'No stagnation point found!'
  STOP
end if

END SUBROUTINE stagnation_point
