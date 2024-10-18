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
! 20.02.2018, TU Dresden, Florian Dexl, WiMi

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

! file module_modify_wake.f90

! This module contains modify_wake subroutine which allows the modification
! of the wake panel orientation to fit the freestream-direction

module module_modify_wake

use module_kind_and_konst

contains

! SUBROUTINE MODIFY_WAKE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine modify_wake( interactive                                     ,&
                          & panel_num                                       ,&
                          & panel_num_no_wake                               ,&
                          & node_num                                        ,&
                          & node1,node2,node3,node4                         ,&
                          & alfa,beta                                       ,&
                          & x_orig, y_orig, z_orig                          ,&
                          & x,y,z                                           )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                       intent(in)                     :: panel_num                  ,&
                                                                   & panel_num_no_wake          ,&
                                                                   & node_num                   ,&
                                                                   & interactive
    integer, dimension(panel_num), intent(in)                     :: node1,node2,node3,node4
    real(kind=kind_float),                          intent(in)    :: alfa,beta                  ,&
                                                                   & x_orig, y_orig, z_orig
    real(kind=kind_float),    dimension(node_num),  intent(inout) :: x,y,z
    ! PRIVATE ======================================================================
    integer                                                       :: i
    logical, dimension(node_num)                                  :: adapted
    real(kind=kind_float)                                         :: dir_x, dir_y, dir_z        ,&
                                                                   & diff_x, diff_y, diff_z     ,&
                                                                   & fac
    ! ==============================================================================
    
    adapted(:) = .false.

    ! Direction vector of flow
    dir_x = cos(alfa)*cos(beta)
    dir_y = sin(beta)
    dir_z = cos(beta)*sin(alfa)
    
    ! Loop over wake panels
    do i=panel_num_no_wake+1,panel_num
        ! Assume, that node 1 and 2 of wake panel are at the wing's trailing edge
        ! and node 3 and 4 are in the Trefftz-plane

        if (adapted(node3(i)) .ne. .true.) then
            ! Projection of node 3 in Trefftz-plane
            diff_x  = x_orig - x(node2(i)); diff_y  = y_orig - y(node2(i)); diff_z  = z_orig - z(node2(i))
            fac     = (dir_x*diff_x + dir_y*diff_y + dir_z*diff_z)/(dir_x**2.d0 + dir_y**2.d0 + dir_z**2.d0)
            x(node3(i)) = x(node2(i)) + dir_x*fac
            y(node3(i)) = y(node2(i)) + dir_y*fac
            z(node3(i)) = z(node2(i)) + dir_z*fac
        
            adapted(node3(i)) = .true.
        end if
        
        if (adapted(node4(i)) .ne. .true.) then
            ! Projection of node 4 in Trefftz-plane
            diff_x  = x_orig - x(node1(i)); diff_y  = y_orig - y(node1(i)); diff_z  = z_orig - z(node1(i))
            fac     = (dir_x*diff_x + dir_y*diff_y + dir_z*diff_z)/(dir_x**2.d0 + dir_y**2.d0 + dir_z**2.d0)
            x(node4(i)) = x(node1(i)) + dir_x*fac
            y(node4(i)) = y(node1(i)) + dir_y*fac
            z(node4(i)) = z(node1(i)) + dir_z*fac
        
            adapted(node4(i)) = .true.
        end if

    end do
    
999 continue
    
    return
    end subroutine modify_wake
    
end module module_modify_wake
