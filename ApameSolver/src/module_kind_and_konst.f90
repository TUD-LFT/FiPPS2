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
! 24.01.2018, TU Dresden, Florian Dexl, WiMi

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

! file module_kind_and_konst.f90

! This module defines and makes available type-parameters and constants.

module module_kind_and_konst

! TYPE PARAMETERS===============================================================
integer, parameter                              :: kind_float_h = 8
integer, parameter                              :: kind_float_l = 4

integer, parameter                              :: kind_float = kind_float_h

integer, parameter                              :: kind_int = 4
! ==============================================================================

! GLOBAL CONSTANTS =============================================================
real(kind=kind_float), parameter                :: pi=4._kind_float*atan(1._kind_float)
! ==============================================================================

! TYPE DEFINITIONS =============================================================
TYPE :: results_type
      integer :: coef       ! coefficients
      integer :: forc       ! forces
      integer :: geom       ! geometry
      integer :: velo       ! velocities
      integer :: pres       ! pressure coefficients
      integer :: cent       ! center points
      integer :: doub       ! doublet values
      integer :: sorc       ! source values
      integer :: velc       ! velocity components
      integer :: mesh       ! mesh characteristics
      integer :: stat       ! static pressure
      integer :: dyna       ! dynamic pressure
      integer :: mano       ! manometer pressure
      integer :: cl_strip   ! stripwise lift coefficient
END TYPE results_type
! ==============================================================================

end module module_kind_and_konst
