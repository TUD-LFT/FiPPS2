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
MODULE functions

use konst_var

implicit none

save

contains

double precision FUNCTION deg_to_rad(deg)
    !
    implicit none

    double precision                                :: deg

    deg_to_rad = deg*pi/180.d0
END FUNCTION deg_to_rad


double precision FUNCTION rad_to_deg(rad)
    !
    implicit none

    double precision                                :: rad

    rad_to_deg = rad*180.d0/pi
END FUNCTION rad_to_deg

FUNCTION cross(v_1, v_2)
    !
    implicit none
    
    double precision, dimension(3)              :: v_1, v_2, cross
    
    cross(1) = v_1(2)*v_2(3) - v_1(3)*v_2(2)
    cross(2) = v_1(3)*v_2(1) - v_1(1)*v_2(3)
    cross(3) = v_1(1)*v_2(2) - v_1(2)*v_2(1)
    
END FUNCTION cross

double precision FUNCTION abs_vector2(vector)
    !
    implicit none
    
    double precision, dimension(2)              :: vector
    
    abs_vector2 = sqrt(vector(1)**2.d0 + vector(2)**2.d0)
    
END FUNCTION abs_vector2

double precision FUNCTION abs_vector(vector)
    !
    implicit none
    
    double precision, dimension(3)              :: vector
    
    abs_vector = sqrt(vector(1)**2.d0 + vector(2)**2.d0 + vector(3)**2.d0)
    
END FUNCTION abs_vector
    
END MODULE functions