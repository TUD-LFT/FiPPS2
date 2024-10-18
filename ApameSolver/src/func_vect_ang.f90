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
! Copyright (C) 2011  Daniel Filkovic

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

! file func_vect_ang.f90

! This function calculates angle between vectors A and B
! Angles must be normalized!!!

subroutine func_vect_ang(A, B, angle)

use module_kind_and_konst

implicit none

! INTENTS IN ===================================================================
real(kind=kind_float), dimension(3), intent(in)  :: A,B
! INTENTS OUT ==================================================================
real(kind=kind_float),               intent(out) :: angle
! PRIVATE ======================================================================
real(kind=kind_float)                            :: cosine
! ==============================================================================

cosine = A(1)*B(1) + A(2)*B(2) + A(3)*B(3);
! check if due to numerical errors we got out of unit circle bounds
if (cosine .ge. 1.0_kind_float) then
    angle = 0.0
elseif (cosine .le. -1.0_kind_float) then
    angle = 3.1415926535897932384626433832795_kind_float
else
    angle = acos(cosine)
endif

return
end subroutine func_vect_ang
