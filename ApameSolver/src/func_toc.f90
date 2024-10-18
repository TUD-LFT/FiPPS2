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
! Copyright (C) 2010-2011  Daniel Filkovic

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

! file func_toc.f90

! This function outputs given text string along with calculated time from
! given start_time up to now.

subroutine func_toc(interactive, start_time, string)

use module_kind_and_konst

implicit none

! INTENTS IN ===================================================================
integer,          intent(in)                  :: interactive
real(kind=kind_float),             intent(in) :: start_time
character(len=*), intent(in)                  :: string
! PRIVATE ======================================================================
integer                                       :: hour, minute
real(kind=kind_float)                         :: total_time, second
character(len=10)                             :: time
! ==============================================================================

! getting current time
call date_and_time(TIME=time)

! reading hour minute and second from "time"
read (time,100) hour
read (time,101) minute
read (time,102) second

! calculating time in seconds and substituting start time
total_time = (hour*3600+minute*60+second) - start_time

! printing text string on screen and/or log file
if (interactive .eq. 1) then
    write (6,103) string,total_time,"s"
    write (2,103) string,total_time,"s"
else
    write (2,103) string,total_time,"s"
endif

! FORMATS ======================================================================
100 format(i2)
101 format(2x,i2)
102 format(4x,F6.3)
103 format(1x,a,F8.3,1x,a)
! ==============================================================================

return
end subroutine func_toc
