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
MODULE netz_variablen

implicit none

save

! Nodal values
integer                                           :: n_nodes
double precision, dimension(:,:), allocatable     :: nodes

! Panel values
integer                                           :: n_panels
double precision, dimension(:), allocatable       :: alpha_i
double precision, dimension(:), allocatable       :: lengths
double precision, dimension(:,:), allocatable     :: coll
double precision, dimension(:,:), allocatable     :: n_vec
double precision, dimension(:,:), allocatable     :: t_vec

!! Wake panel values
!double precision                                  :: alpha_wake
!double precision, dimension(2,2)                  :: wake_nodes

END MODULE netz_variablen
