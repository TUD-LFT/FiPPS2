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
SUBROUTINE unify_airfoil(af_scale, af_rotation)

use konst_var
use netz_variablen

implicit none

! Declaration of variables
double precision, intent(out)                     :: af_scale, af_rotation
double precision                                  :: xmin, xmax
double precision                                  :: rot
double precision, dimension(2)                    :: node_temp
double precision, dimension(2,2)                  :: rot_matrix

integer                                           :: ii
integer                                           :: istart, iend

! Calculation

! Get minimum and maximum of x-values
xmin = huge(1.d0)
xmax = -1.d0*huge(1.d0)
istart = 0
iend   = 0
do ii = 1, n_nodes
    if (nodes(ii,1) .gt. xmax) then 
        xmax = nodes(ii,1)
        iend = ii
    end if
    if (nodes(ii,1) .lt. xmin) then
        xmin = nodes(ii,1)
        istart = ii
    end if
end do

! Set leading edge on (0,0)
node_temp(:) = nodes(istart,:)
do ii = 1, n_nodes
    nodes(ii,:) = nodes(ii,:) - node_temp(:)
end do

! Rotate and scale, so that trailing edge is at (1,0)
af_rotation = atan2(nodes(iend,2) - nodes(istart,2), nodes(iend,1) - nodes(istart,1))
rot_matrix(1,1) =  cos(af_rotation); rot_matrix(1,2) = sin(af_rotation)
rot_matrix(2,1) = -sin(af_rotation); rot_matrix(2,2) = cos(af_rotation)
do ii = 1, n_nodes
    nodes(ii,:) = matmul(rot_matrix(:,:), nodes(ii,:))
end do
! Scaling
af_scale = nodes(iend,1) - nodes(istart,1)
nodes(:,:) = nodes(:,:)/af_scale

END SUBROUTINE unify_airfoil
