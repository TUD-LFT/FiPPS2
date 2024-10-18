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
SUBROUTINE panel_geometry()

use konst_var
use netz_variablen
use functions

implicit none

! Variablendeklaration
integer                                           :: ii
double precision                                  :: max_deltaAlpha_i
!double precision                                  :: wake_length = 40.d0

! Calculation of geometric values
do ii = 1, n_panels
    ! Angle between panel and x-axis
    alpha_i(ii) = -1.d0*atan2((nodes(ii+1,2)-nodes(ii,2)),(nodes(ii+1,1)-nodes(ii,1)))
    if (alpha_i(ii) .lt. 0.d0) then
        alpha_i(ii) = 2.d0*pi+alpha_i(ii)
    end if
    ! Collocation point
    coll(ii,:) = (nodes(ii+1,:) + nodes(ii,:))*0.5d0
    ! Normal vector on panel
    n_vec(ii,1) = sin(alpha_i(ii))
    n_vec(ii,2) = cos(alpha_i(ii))
    ! Tangential vector to panel
    t_vec(ii,1) = cos(alpha_i(ii))
    t_vec(ii,2) = -1.d0*sin(alpha_i(ii))
    ! Length of panel
    lengths(ii) = abs_vector2(nodes(ii+1,:)-nodes(ii,:))
end do

max_deltaAlpha_i = -1.d0*HUGE(1.d0)
do ii = 2, n_panels
  if ((abs(alpha_i(ii) - alpha_i(ii-1)) .gt. max_deltaAlpha_i) .and. (abs(alpha_i(ii) - alpha_i(ii-1)) .lt. pi)) then
    max_deltaAlpha_i = abs(alpha_i(ii) - alpha_i(ii-1))
  end if
end do

write(*,*) 'Maximum panel angle difference: ', rad_to_deg(max_deltaAlpha_i), ' degrees.'

!! Calculation of wake panel geometry
!alpha_wake = (alpha_i(n_panels) - (pi - alpha_i(1)))/2.d0
!wake_nodes(1,:) = (nodes(1,:) + nodes(n_nodes,:))/2.d0
!wake_nodes(2,1) = wake_nodes(1,1) + cos(alpha_wake)*wake_length
!wake_nodes(2,2) = wake_nodes(1,2) - sin(alpha_wake)*wake_length
!write(*,*) 'alpha_1:      ', alpha_i(1)*180.d0/pi
!write(*,*) 'alpha_n:      ', alpha_i(n_panels)*180.d0/pi
!write(*,*) 'alpha_wake:      ', alpha_i(n_panels+1)*180.d0/pi
!write(*,*) 'wake_nodes(2,:): ', wake_nodes(2,:)

END SUBROUTINE panel_geometry
