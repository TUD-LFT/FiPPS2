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
MODULE xfw_konst_var

implicit none

save

! Konstanten
double precision, parameter         :: pi = 3.14159265358979323846d0

! Variablen
logical                             :: viscous
logical                             :: invert_points
logical                             :: cl_given
logical                             :: delta_aoa
double precision                    :: aoa
double precision                    :: cl_target
double precision                    :: reynolds
double precision                    :: Ma
double precision                    :: q_infty
double precision                    :: chordscale
character(len=256)                  :: res_file,res_prefix

! select transition model
integer                             :: n_crit      = 9
double precision                    :: xtr_up      = 1.d0
double precision                    :: xtr_lo      = 1.d0

! Path to Xfoil
character(len=256)                  :: start_xfoil

END MODULE xfw_konst_var
