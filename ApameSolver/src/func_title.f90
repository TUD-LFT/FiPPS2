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

! file func_title.f90

! This function prints program title and version with text for input file request.

subroutine func_title( version                                          ,&
                     & subversion                                       ,&
                     & build_date                                       )
                     
! INTENTS IN ===================================================================
character(len=*),   intent(in)  :: version,subversion
character(len=*),   intent(in)  :: build_date
! PRIVATE ======================================================================
character(len=200)              :: version_line
! ==============================================================================

! assemble version line
version_line = "=                 Version: "//          &
             & trim(version)//                          &
             & "."//                                    &
             & trim(subversion)//                       &
             & "."//                                    &
             & build_date//                             &
             & "                  ="

write(6,"(1x,a)") ""
write(6,"(1x,a)") ""
write(6,"(1x,a)") "========================================================"
write(6,"(1x,a)") "=                                                      ="
write(6,"(1x,a)") "=            APAME - Aircraft Panel Method             ="
write(6,"(1x,a)") "=______________________________________________________="
write(6,"(1x,a)") "=                                                      ="
write(6,"(1x,a)") "=               3D potential flow solver               ="
write(6,"(1x,a)") "=                                                      ="
write(6,"(1x,a)") trim(version_line)
write(6,"(1x,a)") "=                                                      ="
write(6,"(1x,a)") "========================================================"
write(6,"(1x,a)") ""
write(6,"(1x,a)",advance="no") "Please enter input file without .INP extension: "

return
end subroutine func_title

