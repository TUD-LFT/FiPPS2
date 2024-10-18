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
SUBROUTINE drag_squire_young(th_up,H_up,uete_up,th_lo,H_lo,uete_lo,c_drag)

use konst_var

implicit none

! Declaration of variables
double precision, intent(in)                      :: th_up,H_up,uete_up,th_lo,H_lo,uete_lo
double precision, intent(out)                     :: c_drag

! Calculation
c_drag =          2.d0*th_up*uete_up**((min(H_up,2.5d0) + 5.d0)/2.d0)
c_drag = c_drag + 2.d0*th_lo*uete_lo**((min(H_lo,2.5d0) + 5.d0)/2.d0)

END SUBROUTINE drag_squire_young
