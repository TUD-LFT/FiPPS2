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

! file subr_write_coef.f90

! This subroutine prints on screen and/or log file calculated coefficients
! depending on integer interactive.

subroutine subr_write_coef( interactive                                     ,&
                          & case_num                                        ,&
                          & cx,cy,cz                                        ,&
                          & cl,cm,cn                                        ,&
                          & coef_drag,coef_idrag,coef_side,coef_lift        )

use module_kind_and_konst                          

implicit none

! INTENTS IN ===================================================================
integer,                      intent(in)                  :: interactive          ,&
                                                           & case_num
real(kind=kind_float),    dimension(case_num), intent(in) :: cx,cy,cz             ,&
                                                           & cl,cm,cn             ,&
                                                           & coef_drag,coef_idrag ,&
                                                           & coef_side,coef_lift
! PRIVATE ======================================================================
integer                                                   :: i
! ==============================================================================

call func_message( interactive                                              ,&
                 & "    Coefficients in body reference frame:")
call func_new_line(interactive)
if (interactive .eq. 1) then
    call func_message( interactive                                          ,&
                     & "             CX      CY      CZ      CL      CM      CN")
    do i=1,case_num
        write (6,101,advance="no") "        "
        write (6,100,advance="no") cx(i)
        write (6,100,advance="no") cy(i)
        write (6,100,advance="no") cz(i)
        write (6,100,advance="no") cl(i)
        write (6,100,advance="no") cm(i)
        write (6,100,advance="no") cn(i)
        write (2,101,advance="no") "        "
        write (2,100,advance="no") cx(i)
        write (2,100,advance="no") cy(i)
        write (2,100,advance="no") cz(i)
        write (2,100,advance="no") cl(i)
        write (2,100,advance="no") cm(i)
        write (2,100,advance="no") cn(i)
    enddo
else
    call func_message( interactive                                          ,&
                     & "             CX      CY      CZ      CL      CM      CN")
    do i=1,case_num
        write (2,101,advance="no") "        "
        write (2,100,advance="no") cx(i)
        write (2,100,advance="no") cy(i)
        write (2,100,advance="no") cz(i)
        write (2,100,advance="no") cl(i)
        write (2,100,advance="no") cm(i)
        write (2,100,advance="no") cn(i)
    enddo
endif
call func_new_line(interactive)
call func_new_line(interactive)

call func_message( interactive                                              ,&
                 & "    Coefficients in aerodynamic reference frame:")
call func_new_line(interactive)
if (interactive .eq. 1) then
    call func_message( interactive                                          ,&
                     & "             CD     CDi      CK      CL")
    do i=1,case_num
        write (6,101,advance="no") "        "
        write (6,100,advance="no") coef_drag(i)
        write (6,100,advance="no") coef_idrag(i)
        write (6,100,advance="no") coef_side(i)
        write (6,100,advance="no") coef_lift(i)
        write (2,101,advance="no") "        "
        write (2,100,advance="no") coef_drag(i)
        write (2,100,advance="no") coef_idrag(i)
        write (2,100,advance="no") coef_side(i)
        write (2,100,advance="no") coef_lift(i)
    enddo
else
    call func_message( interactive                                          ,&
                     & "             CD     CDi      CK      CL")
    do i=1,case_num
        write (2,101,advance="no") "        "
        write (2,100,advance="no") coef_drag(i)
        write (2,100,advance="no") coef_idrag(i)
        write (2,100,advance="no") coef_side(i)
        write (2,100,advance="no") coef_lift(i)
    enddo
endif

! FORMATS ======================================================================
100 format(1x,F7.4)
101 format(/,a)
! ==============================================================================

return
end subroutine subr_write_coef
