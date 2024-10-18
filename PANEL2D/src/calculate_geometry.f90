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
SUBROUTINE calculate_geometry()

use konst_var
use netz_variablen
use functions

implicit none

! Variablendeklaration
double precision, dimension(:), allocatable       :: ss
integer                                           :: ii
!integer                                           :: n_nodes_new
!double precision, dimension(:), allocatable       :: s
!double precision, dimension(:,:), allocatable     :: nodes_temp

allocate(ss(n_nodes))

! Running length along airfoil
ss(1) = 0.d0
do ii = 2, n_nodes
    ss(ii) = ss(ii-1) + abs_vector2(nodes(ii,:) - nodes(ii-1,:))
end do

! ! Verfeinerung durch Spline-Interpolation
! n_nodes_new = 400
! allocate(nodes_temp(n_nodes_new,2),s(n_nodes_new))
! do ii = 1, n_nodes_new
!    s(ii) = dble(ii-1)/(n_nodes_new-1)*ss(n_nodes)
! end do
! call splineInterpol(n_nodes,n_nodes_new,ss(:),nodes(:,1),s(:),nodes_temp(:,1))
! call splineInterpol(n_nodes,n_nodes_new,ss(:),nodes(:,2),s(:),nodes_temp(:,2))
! deallocate(nodes)
! allocate(nodes(n_nodes_new,2))
! nodes(:,:) = nodes_temp(:,:)
! deallocate(nodes_temp)
! n_nodes = n_nodes_new
! n_panels = n_nodes-1

! Allocate arrays for geometric sizes
allocate( n_vec(n_panels,2),     &
        & t_vec(n_panels,2),     &
        & lengths(n_panels),     &
        & alpha_i(n_panels),     &
        & coll(n_panels,2)       )

END SUBROUTINE calculate_geometry
