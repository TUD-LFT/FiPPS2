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

! file subr_write_ptot.f90

! This subroutine prints on screen and/or log file total number of
! panels depending on integer interactive.

subroutine subr_write_ptot( interactive                                     ,&
                          & panel_num                                       ,&
                          & panel_num_no_wake                               )

implicit none

! INTENTS IN ===================================================================
integer, intent(in) :: interactive                                          ,&
                     & panel_num                                            ,&
                     & panel_num_no_wake
! ==============================================================================

call func_message( interactive                                              ,&
                 & "    Total number of panels:    ")
if (interactive .eq. 1) then
    write (6,100) panel_num
    write (2,100) panel_num
else
    write (2,100) panel_num
endif
call func_message( interactive                                              ,&
                 & "    ...excluding wake panels:  ")
if (interactive .eq. 1) then
    write (6,100) panel_num_no_wake
    write (2,100) panel_num_no_wake
else
    write (2,100) panel_num_no_wake
endif
call func_new_line(interactive)

! FORMATS ======================================================================
100 format(i7)
! ==============================================================================

return
end subroutine subr_write_ptot
