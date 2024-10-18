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
MODULE xfw_functions

use xfw_konst_var

implicit none

save

contains

double precision FUNCTION abs_vector2(vector)
    !
    implicit none
    
    double precision, dimension(2)              :: vector
    
    abs_vector2 = sqrt(vector(1)**2.d0 + vector(2)**2.d0)
    
END FUNCTION abs_vector2
    
END MODULE xfw_functions
