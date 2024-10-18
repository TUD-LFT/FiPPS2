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

! file subr_write_case.f90

! This subroutine prints on screen and/or log file requested angles of attack
! (alfa) and sideslip angles (beta) in degrees depending on integer interactive.

subroutine subr_write_case(interactive, case_num, alfa, beta, trim_cl, target_cl)

use module_kind_and_konst

implicit none

! INTENTS IN ===================================================================
integer,                      intent(in)                  :: interactive, case_num ,&
                                                           & trim_cl
real(kind=kind_float),    dimension(case_num), intent(in) :: alfa
real(kind=kind_float),    dimension(case_num), intent(in) :: beta
real(kind=kind_float),                         intent(in) :: target_cl
! PRIVATE ======================================================================
integer                                                   :: i
real(kind=kind_float),    dimension(case_num)             :: alfa_deg
real(kind=kind_float),    dimension(case_num)             :: beta_deg
! ==============================================================================

! converting angles to degrees
alfa_deg = alfa/pi*180
beta_deg = beta/pi*180
    
if (trim_cl .eq. 0) then
    ! printing angles of attack
    call func_message( interactive                                              ,&
                    & "    Angles of interest [deg]:")
    call func_new_line(interactive)
    call func_message( interactive                                              ,&
                    & "        Angles of attack:")
    if (interactive .eq. 1) then
        do i=1,case_num
            write (6,100,advance="no") alfa_deg(i)
            write (2,100,advance="no") alfa_deg(i)
        enddo
    else
        do i=1,case_num
            write (2,100,advance="no") alfa_deg(i)
        enddo
    endif
else
    ! printing target lift coefficient
    call func_message( interactive                                              ,&
                    & "    Trimming to target lift coefficient:")
    call func_new_line(interactive)
    call func_message( interactive                                              ,&
                    & "        Target lift coefficient:")
    if (interactive .eq. 1) then
        write (6,100,advance="no") target_cl
    endif
    write (2,100,advance="no") target_cl
end if

! printing sideslip angles
call func_new_line(interactive)
call func_message( interactive                                              ,&
                 & "        Sideslip angles: ")
if (interactive .eq. 1) then
    do i=1,case_num
        write (6,100,advance="no") beta_deg(i)
        write (2,100,advance="no") beta_deg(i)
    enddo
else
    do i=1,case_num
        write (2,100,advance="no") beta_deg(i)
    enddo
endif
call func_new_line(interactive)

! FORMATS ======================================================================
100 format(1x,f5.1)
! ==============================================================================

return
end subroutine subr_write_case
