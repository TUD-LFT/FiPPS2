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

! file func_IDW.f90

! This function calculates weights for Inverse Distance Weighting method.

subroutine weight_factors( weights                                          ,&
                         & node1                                            ,&
                         & node2                                            ,&
                         & node3                                            ,&
                         & node4                                            ,&
                         & x                                                ,&
                         & y                                                ,&
                         & z                                                ,&
                         & cx                                               ,&
                         & cy                                               ,&
                         & cz                                               ,&
                         & node_num                                         ,&
                         & panel_num                                        ,&
                         & panel_type)

use module_kind_and_konst

implicit none

! INTENTS IN ===================================================================
integer,                            intent(in)                       :: node_num,panel_num
integer, dimension(panel_num),      intent(in)                       :: panel_type
integer, dimension(panel_num),      intent(in)                       :: node1,node2,node3,node4
real(kind=kind_float),    dimension(panel_num),      intent(in)      :: cx,cy,cz
real(kind=kind_float),    dimension(node_num),       intent(in)      :: x,y,z
! INTENTS OUT ==================================================================
real(kind=kind_float),    dimension(panel_num,4),    intent(inout)   :: weights
! PRIVATE ======================================================================
integer                                                              :: i
! ==============================================================================

do i=1,panel_num
    ! calculate inverse distance for each node of current element
    weights(i,1)=1/( (x(node1(i)) - cx(i))**2   +&
                   & (y(node1(i)) - cy(i))**2   +&
                   & (z(node1(i)) - cz(i))**2   )
    
    weights(i,2)=1/( (x(node2(i)) - cx(i))**2   +&
                   & (y(node2(i)) - cy(i))**2   +&
                   & (z(node2(i)) - cz(i))**2   )
    
    weights(i,3)=1/( (x(node3(i)) - cx(i))**2   +&
                   & (y(node3(i)) - cy(i))**2   +&
                   & (z(node3(i)) - cz(i))**2   )
    
    if ( panel_type(i) .eq. 1   .or.&
       & panel_type(i) .eq. 10  .or.&
       & panel_type(i) .eq. 20  )then
        weights(i,4)=1/( (x(node4(i)) - cx(i))**2   +&
                       & (y(node4(i)) - cy(i))**2   +&
                       & (z(node4(i)) - cz(i))**2   )
    endif
enddo

return
end subroutine weight_factors
