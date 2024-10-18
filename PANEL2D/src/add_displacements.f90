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
SUBROUTINE add_displacements(displacements, scle)

use konst_var
use netz_variablen

implicit none

! Declaration of variables
double precision, dimension(:,:), intent(in)      :: displacements
double precision, intent(in)                      :: scle

! Calculation

if (size(displacements,1) .ne. (n_nodes - 1)) then
    write(*,*) 'Wrong number of displacements given!'
    stop
end if
if (invert_points .eqv. .false.) then
    nodes(1:(n_nodes-1),:) = nodes(1:(n_nodes-1),:) + displacements(:,1:2)*scle
    nodes(n_nodes,:) = nodes(n_nodes,:) + displacements(1,1:2)*scle
else
    nodes(2:n_nodes,:) = nodes(2:n_nodes,:) + displacements((n_nodes-1):1:-1,1:2)*scle
    nodes(1,:) = nodes(1,:) + displacements(1,1:2)*scle
end if

END SUBROUTINE add_displacements
