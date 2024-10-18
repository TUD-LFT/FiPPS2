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
! 31.01.2018, TU Dresden, Florian Dexl, WiMi

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

! file subr_vtk_results.f90

! This subroutine writes result data to a vtk (Visualization Toolkit) file

subroutine subr_vtk_results( interactive                                    ,&
                       & version                                            ,&
                       & subversion                                         ,&
                       & build_date                                         ,&
                       & file_name                                          ,&
                       & speed                                              ,&
                       & ro                                                 ,&
                       & p_ref                                              ,&
                       & mach                                               ,&
                       & alfa                                               ,&
                       & beta                                               ,&
                       & wing_span                                          ,&
                       & mac                                                ,&
                       & wing_surf                                          ,&
                       & origin                                             ,&
                       & method                                             ,&
                       & error                                              ,&
                       & coplanarity_angle                                  ,&
                       & farfield                                           ,&
                       & collcalc                                           ,&
                       & velorder                                           ,&
                       & case_num                                           ,&
                       & node_num                                           ,&
                       & panel_num                                          ,&
                       & panel_num_no_wake                                  ,&
                       & panel_type                                         ,&
                       & x,y,z                                              ,&
                       & node1,node2,node3,node4                            ,&
                       & elem1,elem2,elem3,elem4                            ,&
                       & cx,cy,cz                                           ,&
                       & cp                                                 ,&
                       & v                                                  ,&
                       & vx,vy,vz                                           ,&
                       & p_dyna                                             ,&
                       & p_mano                                             ,&
                       & p_stat                                             ,&
                       & gama                                               ,&
                       & sigma                                              ,&
                       & S                                                  ,&
                       & FF                                                 ,&
                       & n1,n2,n3                                           ,&
                       & l1,l2,l3                                           ,&
                       & p1,p2,p3                                           ,&
                       & coef_x,coef_y,coef_z                               ,&
                       & coef_l,coef_m,coef_n                               ,&
                       & coef_drag,coef_idrag,coef_side,coef_lift           ,&
                       & cl_strip, cllength_strip, elem_strip_id            ,&
                       & Fx,Fy,Fz                                           ,&
                       & Fl,Fm,Fn                                           ,&
                       & Fdrag,Fidrag,Fside,Flift                           ,&
                       & results                                            )

use module_kind_and_konst

implicit none
! 
! INTENTS IN ===================================================================
integer,                                   intent(in)                  :: interactive        ,&
                                                                        & case_num           ,&
                                                                        & collcalc           ,&
                                                                        & velorder           ,&
                                                                        & method             ,&
                                                                        & node_num           ,&
                                                                        & panel_num          ,&
                                                                        & panel_num_no_wake
integer,    dimension(panel_num),          intent(in)                  :: panel_type         ,&
                                                                        & node1,node2        ,&
                                                                        & node3,node4        ,&
                                                                        & elem1,elem2        ,&
                                                                        & elem3,elem4
real(kind=kind_float),                                      intent(in) :: speed              ,&
                                                                        & ro                 ,&
                                                                        & p_ref              ,&
                                                                        & mach               ,&
                                                                        & wing_span          ,&
                                                                        & mac                ,&
                                                                        & wing_surf          ,&
                                                                        & error              ,&
                                                                        & coplanarity_angle  ,&
                                                                        & farfield           ,&
                                                                        & origin(3)
real(kind=kind_float),       dimension(case_num),           intent(in) :: alfa,beta          ,&
                                                                        & Fx,Fy,Fz           ,&
                                                                        & Fl,Fm,Fn           ,&
                                                                        & Fdrag,Fidrag       ,&
                                                                        & Fside,Flift        ,&
                                                                        & coef_x,coef_y,coef_z,&
                                                                        & coef_l,coef_m,coef_n,&
                                                                        & coef_drag          ,&
                                                                        & coef_idrag         ,&
                                                                        & coef_side          ,&
                                                                        & coef_lift
real(kind=kind_float),       dimension(node_num),           intent(in) :: x,y,z
real(kind=kind_float),       dimension(panel_num),          intent(in) :: S                  ,&
                                                                        & FF                 ,&
                                                                        & cx,cy,cz           ,&
                                                                        & n1,n2,n3           ,&
                                                                        & l1,l2,l3           ,&
                                                                        & p1,p2,p3
real(kind=kind_float),       dimension(panel_num_no_wake,case_num), intent(in) :: cp         ,&
                                                                        & v                  ,&
                                                                        & vx,vy,vz           ,&
                                                                        & gama               ,&
                                                                        & sigma              ,&
                                                                        & p_dyna             ,&
                                                                        & p_mano             ,&
                                                                        & p_stat
real(kind=kind_float),       dimension(panel_num-panel_num_no_wake,case_num), intent(in) :: cl_strip, cllength_strip
integer              ,       dimension(panel_num_no_wake),  intent(in) :: elem_strip_id
character(len=100),                        intent(in)                  :: file_name
character(len=3),                          intent(in)                  :: version           ,&
                                                                        & subversion
character(len=6),                          intent(in)                  :: build_date
TYPE(results_type)                                                     :: results
! PRIVATE ======================================================================
integer                                                                :: i,j                ,&
                                                                        & year               ,&
                                                                        & month              ,&
                                                                        & day                ,&
                                                                        & hour               ,&
                                                                        & minute
integer                                                                :: n_tri,n_quad
integer                                                                :: unit_vtk
character(len=8)                                                       :: date
character(len=10)                                                      :: time
character(len=10)                                                      :: formatter

parameter (unit_vtk=30)
! ==============================================================================

! open results vtk file for writing
open(unit_vtk,file = trim(file_name)//".vtk",status="replace")

! write header
write(unit_vtk,'(A26)') '# vtk DataFile Version 4.2'
write(formatter,'("(A",I0,")")') len_trim(file_name)
write(unit_vtk,formatter) trim(file_name)
write(unit_vtk,'(A5)') 'ASCII'
write(unit_vtk,'(A25)') 'DATASET UNSTRUCTURED_GRID'

! write nodes
write(unit_vtk,'(A6,X,I0,X,A6)') 'POINTS', node_num, 'DOUBLE'
do i=1,node_num
    write(unit_vtk,'(3E24.16)') x(i), y(i), z(i)
end do

! count quadrilateral panels (incl. wake panels)
n_quad = 0
do i=1,panel_num
    if (panel_type(i) .eq. 1 .or. panel_type(i) .eq. 10 .or. panel_type(i) .eq. 20) then
        n_quad = n_quad+1
    end if
end do

! count triangular panels (incl. wake panels)
n_tri = 0
do i=1,panel_num
    if (panel_type(i) .eq. 2 .or. panel_type(i) .eq. 11 .or. panel_type(i) .eq. 21) then
        n_tri = n_tri+1
    end if
end do

! writing panels
write(unit_vtk,'(A5,X,I0,X,I0)') 'CELLS', panel_num, (n_tri*4 + n_quad*5)
do i=1,panel_num
    if (panel_type(i) .eq. 1 .or. panel_type(i) .eq. 10 .or. panel_type(i) .eq. 20) then
        write(unit_vtk,'(I0,X,I0,X,I0,X,I0,X,I0)') 4,(node1(i)-1),(node2(i)-1),(node3(i)-1),(node4(i)-1)
    else if (panel_type(i) .eq. 2 .or. panel_type(i) .eq. 11 .or. panel_type(i) .eq. 21) then
        write(unit_vtk,'(I0,X,I0,X,I0,X,I0)') 3,(node1(i)-1),(node2(i)-1),(node3(i)-1)
    end if
end do

! writing cell types
write(unit_vtk,'(A10,X,I0)') 'CELL_TYPES', panel_num
do i=1,panel_num
    if (panel_type(i) .eq. 1 .or. panel_type(i) .eq. 10 .or. panel_type(i) .eq. 20) then
        write(unit_vtk,'(I0)') 9
    else if (panel_type(i) .eq. 2 .or. panel_type(i) .eq. 11 .or. panel_type(i) .eq. 21) then
        write(unit_vtk,'(I0)') 5
    end if
end do

! writing panel types
write(unit_vtk,'(A9,X,I0)') 'CELL_DATA', panel_num
write(unit_vtk,'(A7,X,A4,X,A5)') 'SCALARS', 'TYPE', 'INT 1'
write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
do i=1,panel_num
    write(unit_vtk,'(I0)') panel_type(i)
end do

! writing pressure coefficients on panels
if (results%pres .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A8,I0,X,A8)') 'SCALARS', 'CP_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') cp(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing velocities on panels
if (results%velo .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A9,I0,X,A8)') 'SCALARS', 'VEL_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') v(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing dipole strengths on panels
if (results%doub .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A12,I0,X,A8)') 'SCALARS', 'DIPOLE_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') gama(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') gama(elem1(j),i) - gama(elem2(j),i)
        end do
    end do
end if

! writing source strengths on panels
! if requested and if available (doublet/source combination only)
if (results%sorc .eq. 1 .and. method .ne. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A12,I0,X,A8)') 'SCALARS', 'SOURCE_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') sigma(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing static pressure on panels
if (results%stat .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A14,I0,X,A8)') 'SCALARS', 'P_STATIC_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') p_stat(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing dynamic pressure on panels
if (results%dyna .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A15,I0,X,A8)') 'SCALARS', 'P_DYNAMIC_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') p_dyna(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing manometer pressure on panels
if (results%mano .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A23,I0,X,A8)') 'SCALARS', 'P_MANOMETER_CASE_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            write(unit_vtk,'(E24.16)') p_mano(j,i)
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing stripwise lift coefficients on panels
if (results%cl_strip .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A9,I0,X,A8)') 'SCALARS', 'CL_STRIP_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            if (elem_strip_id(j) .ne. 0) then
                write(unit_vtk,'(E24.16)') cl_strip(elem_strip_id(j),i)
            else
                write(unit_vtk,'(E24.16)') 0._kind_float
            end if
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing stripwise lift coefficients multiplied with local chords on panels
if (results%cl_strip .eq. 1) then
    do i=1,case_num
        write(unit_vtk,'(A7,X,A19,I0,X,A8)') 'SCALARS', 'CL_STRIP*LOC_CHORD_', i, 'DOUBLE 1'
        write(unit_vtk,'(A12,X,A7)') 'LOOKUP_TABLE', 'DEFAULT'
        do j=1,panel_num_no_wake
            if (elem_strip_id(j) .ne. 0) then
                write(unit_vtk,'(E24.16)') cllength_strip(elem_strip_id(j),i)
            else
                write(unit_vtk,'(E24.16)') 0._kind_float
            end if
        end do
        do j=panel_num_no_wake+1,panel_num
            write(unit_vtk,'(E24.16)') 0._kind_float
        end do
    end do
end if

! writing mesh characteristics
if (results%mesh .eq. 1) then
    write(unit_vtk,'(A7,X,A6,X,A6)') 'NORMALS', 'NORMAL', 'DOUBLE'
    do i=1,panel_num
        write(unit_vtk,'(3E24.16)') n1(i), n2(i), n3(i)
    end do

    write(unit_vtk,'(A7,X,A8,X,A6)') 'VECTORS', 'L_VECTOR', 'DOUBLE'
    do i=1,panel_num
        write(unit_vtk,'(3E24.16)') l1(i), l2(i), l3(i)
    end do

    write(unit_vtk,'(A7,X,A8,X,A6)') 'VECTORS', 'P_VECTOR', 'DOUBLE'
    do i=1,panel_num
        write(unit_vtk,'(3E24.16)') p1(i), p2(i), p3(i)
    end do
end if

! close results vtk file
close(unit_vtk)

return
end subroutine subr_vtk_results
