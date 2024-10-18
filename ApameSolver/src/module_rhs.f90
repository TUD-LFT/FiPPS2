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

! file module_rhs.f90

! This function calculates right-hand side vector using grid data and free stream
! velocity components. If combination of constant source/doublet is used then first
! source values are calculated and then this vector is multiplied with matrix b.
! If constant doublet method is used then right-hand side represents free stream
! velocity potential

module module_rhs

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                           :: rhs_err
real(kind=kind_float), allocatable, dimension(:,:) :: sigma, rhs
! ==============================================================================

contains

    subroutine rhs_calc( b                                                  ,&
                       & method                                             ,&
                       & panel_num                                          ,&
                       & panel_num_no_wake                                  ,&
                       & interactive                                        ,&
                       & case_num                                           ,&
                       & speed_x,speed_y,speed_z                            ,&
                       & n1,n2,n3                                           ,&
                       & cx,cy,cz                                           )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                                intent(in)                  :: method            ,&
                                                                         & panel_num         ,&
                                                                         & panel_num_no_wake ,&
                                                                         & interactive       ,&
                                                                         & case_num
    real(kind=kind_float),    dimension(panel_num),          intent(in) :: n1,n2,n3          ,&
                                                                         & cx,cy,cz
    real(kind=kind_float),    dimension(case_num),           intent(in) :: speed_x,speed_y,speed_z
    real(kind=kind_float),    dimension( panel_num_no_wake                                   ,&
                                       & panel_num_no_wake), intent(in) :: b
    ! PRIVATE ======================================================================
    integer                                                             :: i,j               ,&
                                                                         & alloc_stat
    ! ==============================================================================
    
    rhs_err=0
    
    ! if constant source/doublet method defined
    if (method .eq. 0) then
        ! calculating source strengths
        allocate ( sigma(panel_num_no_wake,case_num)    ,&
                 & rhs(panel_num_no_wake,case_num)      ,&
                 & stat = alloc_stat                    )
        if (alloc_stat .ne. 0) then
            call func_message( interactive                                  ,&
                 & "    ERROR: Not enough memory for allocating source and RHS vectors")
            rhs_err = 1
            goto 999
        endif
        
        do i=1,panel_num_no_wake
            do j=1,case_num
                sigma(i,j) = n1(i)*speed_x(j) + n2(i)*speed_y(j) + n3(i)*speed_z(j)
            enddo
        enddo
        
        ! calculating "right hand side" by multiplying source influence matrix (b) with
        ! source strengths (sigma)
        if (kind_float .eq. kind_float_l) then
            call sgemm( "N"                                                     ,&
                      & "N"                                                     ,&
                      & panel_num_no_wake                                       ,&
                      & case_num                                                ,&
                      & panel_num_no_wake                                       ,&
                      & 1.0e0_kind_float                                        ,&
                      & b                                                       ,&
                      & panel_num_no_wake                                       ,&
                      & sigma                                                   ,&
                      & panel_num_no_wake                                       ,&
                      & 0.0e0_kind_float                                        ,&
                      & rhs                                                     ,&
                      & panel_num_no_wake                                       )
        else
            call dgemm( "N"                                                     ,&
                      & "N"                                                     ,&
                      & panel_num_no_wake                                       ,&
                      & case_num                                                ,&
                      & panel_num_no_wake                                       ,&
                      & 1.0e0_kind_float                                        ,&
                      & b                                                       ,&
                      & panel_num_no_wake                                       ,&
                      & sigma                                                   ,&
                      & panel_num_no_wake                                       ,&
                      & 0.0e0_kind_float                                        ,&
                      & rhs                                                     ,&
                      & panel_num_no_wake                                       )
        end if
    
    ! if constant doublet method defined
    elseif (method .eq. 1) then
        ! calculating free stream potential
        allocate ( rhs(panel_num_no_wake,case_num)      ,&
                 & stat = alloc_stat                    )
        if (alloc_stat .ne. 0) then
            call func_message( interactive                                  ,&
                 & "    ERROR: Not enough memory for allocating RHS vector")
            rhs_err = 1
            goto 999
        endif
        
        do i=1,panel_num_no_wake
            do j=1,case_num
                rhs(i,j) = cx(i)*speed_x(j) + cy(i)*speed_y(j) + cz(i)*speed_z(j)
            enddo
        enddo
    endif
    
999 continue
    
    return
    end subroutine rhs_calc

end module module_rhs
