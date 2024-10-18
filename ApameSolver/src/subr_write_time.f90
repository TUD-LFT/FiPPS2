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
! Copyright (C) 2010  Daniel Filkovic

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %                                                      %
! %            APAME - Aircraft Panel Method             %
! %______________________________________________________%
! %                                                      %
! %              3D potential flow solver                %
! %                                                      %
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This file is part of APAME.

! APAME is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! APAME is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with APAME.  If not, see <http://www.gnu.org/licenses/>.

! file subr_write_time.f90

! This subroutine prints on screen and/or log file date and time
! depending on integer interactive.

subroutine subr_write_time(interactive)

implicit none

! INTENTS IN ===================================================================
integer,      intent(in) :: interactive
! PRIVATE ======================================================================
character(8)             :: date
character(10)            :: time
integer                  :: year,month,day,hour,minute
! ==============================================================================

! finding date and time
call date_and_time(DATE=date,TIME=time)
read (date,100) year
read (date,101) month
read (date,102) day
read (time,103) hour
read (time,104) minute

! printing date and time
if (interactive .eq. 1) then
    write (6,106) day,month,year,hour,minute
    write (2,106) day,month,year,hour,minute
else
    write (2,106) day,month,year,hour,minute
endif

! FORMATS ======================================================================
100 format(i4)
101 format(4x,i2)
102 format(6x,i2)
103 format(i2)
104 format(2x,i2)
106 format(1x,"Date:",1x,i2.2,"/",i2.2,"/",i4,2x,"Time:",1x,i2.2,":",i2.2)
! ==============================================================================

return
end subroutine subr_write_time
