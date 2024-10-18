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
SUBROUTINE sort_nodes()

use konst_var
use netz_variablen
use functions

implicit none

! Declaration of variables
double precision                                  :: z_side1, z_side2
double precision, dimension(:,:), allocatable     :: nodes_temp
integer                                           :: le_pos
integer                                           :: ii

! Calculation

! If necessary, invert order of nodes
! For later calculations, the order of nodes must run from the
! trailing edge via the airfoil's lower surface to the leading
! edge and via the airfoil's upper surface back to the trailing
! edge
! Find position of the leading edge
le_pos  = minloc(nodes(:,1), 1)
! Check whether node order has to be inverted
z_side1 = sum(nodes(1:le_pos,2))
z_side2 = sum(nodes(le_pos+1:,2))
if (z_side1 .gt. z_side2) then
    invert_points = .true.
else
    invert_points = .false.
end if

if (invert_points .eqv. .true.) then
    ! Invert node order
    allocate(nodes_temp(n_nodes,2))
    nodes_temp(:,:) = nodes(:,:)
    do ii = 1, n_nodes
        nodes(ii,:) = nodes_temp(n_nodes-ii+1,:)
    end do
    deallocate(nodes_temp)
end if

END SUBROUTINE sort_nodes
