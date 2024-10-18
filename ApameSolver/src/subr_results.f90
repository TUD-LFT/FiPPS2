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
! Copyright (C) 2010-2012  Daniel Filkovic

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

! file subr_results.f90

! This subroutine writes result data to results file

subroutine subr_results( interactive                                        ,&
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
                       & cl_strip, cllength_strip, strip_center             ,&
                       & Fx,Fy,Fz                                           ,&
                       & Fl,Fm,Fn                                           ,&
                       & Fdrag,Fidrag,Fside,Flift                           ,&
                       & results                                            )

use module_kind_and_konst

implicit none

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
real(kind=kind_float),       dimension(panel_num-panel_num_no_wake,3),        intent(in) :: strip_center
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
character(len=8)                                                       :: date
character(len=10)                                                      :: time
! ==============================================================================

! check if source distribution requested, but does not exist
if (method .eq. 1 .and. results%sorc .eq. 1) then
    call func_new_line(interactive)
    call func_message( interactive                                          ,&
                     & " WARNING: SOURCE distribution not available for this method...")
    call func_new_line(interactive)
endif

! open results file for writing
open(unit=3,file = trim(file_name)//".res",status="replace")
call date_and_time(DATE=date,TIME=time)
read (date,100) year
read (date,101) month
read (date,102) day
read (time,103) hour
read (time,104) minute

write (3,206) "APAME results file"
write (3,106) day,month,year,hour,minute
write (3,206) "VERSION"
write (3,206) trim(version)//"."//trim(subversion)//"."//build_date
write (3,206) "PANEL_NUM"
write (3,208) panel_num
write (3,206) "PANEL_NUM_NO_WAKE"
write (3,208) panel_num_no_wake
write (3,206) "NODE_NUM"
write (3,208) node_num
write (3,206) "AIRSPEED"
write (3,209) speed
write (3,206) "DENSITY"
write (3,209) ro
write (3,206) "PRESSURE"
write (3,209) p_ref
write (3,206) "MACH"
write (3,209) mach
write (3,206) "CASE_NUM"
write (3,208) case_num
do i=1,case_num
    write (3,209,advance="no") alfa(i)
enddo
write (3,105) ""
do i=1,case_num
    write (3,209,advance="no") beta(i)
enddo
write (3,105) ""
write (3,206) "WINGSPAN"
write (3,209) wing_span
write (3,206) "MAC"
write (3,209) mac
write (3,206) "SURFACE"
write (3,209) wing_surf
write (3,206) "ORIGIN"
do i=1,3
    write (3,209,advance="no") origin(i)
enddo
write (3,105) ""
write (3,206) "METHOD"
write (3,208) method
write (3,206) "ERROR"
write (3,209) error
write (3,206) "COPL_ANG"
write (3,209) coplanarity_angle
write (3,206) "FARFIELD"
write (3,209) farfield
write (3,206) "COLLCALC"
write (3,208) collcalc
write (3,206) "VELORDER"
write (3,208) velorder

if (results%coef .eq. 1) then
    write (3,206) "CX"
    do i=1,case_num
        write (3,209,advance="no") coef_x(i)
    enddo
    write (3,105) ""
    write (3,206) "CY"
    do i=1,case_num
        write (3,209,advance="no") coef_y(i)
    enddo
    write (3,105) ""
    write (3,206) "CZ"
    do i=1,case_num
        write (3,209,advance="no") coef_z(i)
    enddo
    write (3,105) ""
    write (3,206) "CL"
    do i=1,case_num
        write (3,209,advance="no") coef_l(i)
    enddo
    write (3,105) ""
    write (3,206) "CM"
    do i=1,case_num
        write (3,209,advance="no") coef_m(i)
    enddo
    write (3,105) ""
    write (3,206) "CN"
    do i=1,case_num
        write (3,209,advance="no") coef_n(i)
    enddo
    write (3,105) ""
    write (3,206) "CDRAG"
    do i=1,case_num
        write (3,209,advance="no") coef_drag(i)
    enddo
    write (3,105) ""
    write (3,206) "CDRAGi"
    do i=1,case_num
        write (3,209,advance="no") coef_idrag(i)
    enddo
    write (3,105) ""
    write (3,206) "CSIDE"
    do i=1,case_num
        write (3,209,advance="no") coef_side(i)
    enddo
    write (3,105) ""
    write (3,206) "CLIFT"
    do i=1,case_num
        write (3,209,advance="no") coef_lift(i)
    enddo
    write (3,105) ""
endif

if (results%forc .eq. 1) then
    write (3,206) "FX"
    do i=1,case_num
        write (3,209,advance="no") Fx(i)
    enddo
    write (3,105) ""
    write (3,206) "FY"
    do i=1,case_num
        write (3,209,advance="no") Fy(i)
    enddo
    write (3,105) ""
    write (3,206) "FZ"
    do i=1,case_num
        write (3,209,advance="no") Fz(i)
    enddo
    write (3,105) ""
    write (3,206) "FL"
    do i=1,case_num
        write (3,209,advance="no") Fl(i)
    enddo
    write (3,105) ""
    write (3,206) "FM"
    do i=1,case_num
        write (3,209,advance="no") Fm(i)
    enddo
    write (3,105) ""
    write (3,206) "FN"
    do i=1,case_num
        write (3,209,advance="no") Fn(i)
    enddo
    write (3,105) ""
    write (3,206) "FDRAG"
    do i=1,case_num
        write (3,209,advance="no") Fdrag(i)
    enddo
    write (3,105) ""
    write (3,206) "FDRAGi"
    do i=1,case_num
        write (3,209,advance="no") Fidrag(i)
    enddo
    write (3,105) ""
    write (3,206) "FSIDE"
    do i=1,case_num
        write (3,209,advance="no") Fside(i)
    enddo
    write (3,105) ""
    write (3,206) "FLIFT"
    do i=1,case_num
        write (3,209,advance="no") Flift(i)
    enddo
    write (3,105) ""
endif

! writing geometry (nodes and panels)
if (results%geom .eq. 1) then
    ! nodes
    write (3,206) "NODES"
    do i=1,node_num
        write (3,209,advance="no") x(i)
        write (3,209,advance="no") y(i)
        write (3,209,advance="no") z(i)
        write (3,105) ""
    enddo
    
    ! panels
    write (3,206) "PANELS"
    do i=1,panel_num
        write (3,208,advance="no") panel_type(i)
        write (3,208,advance="no") node1(i)
        write (3,208,advance="no") node2(i)
        write (3,208,advance="no") node3(i)
        write (3,208,advance="no") node4(i)
        write (3,208,advance="no") elem1(i)
        write (3,208,advance="no") elem2(i)
        write (3,208,advance="no") elem3(i)
        write (3,208,advance="no") elem4(i)
        write (3,105) ""
    enddo
endif

! writing mesh characteristics
if (results%mesh .eq. 1) then
    write (3,206) "S"
    do i=1,panel_num
        write (3,209) S(i)
    enddo
    
    write (3,206) "FF"
    do i=1,panel_num
        write (3,209) FF(i)
    enddo
    
    write (3,206) "NORMAL"
    do i=1,panel_num
        write (3,209,advance="no") n1(i)
        write (3,209,advance="no") n2(i)
        write (3,209,advance="no") n3(i)
        write (3,105) ""
    enddo
    
    write (3,206) "L_VECTOR"
    do i=1,panel_num
        write (3,209,advance="no") l1(i)
        write (3,209,advance="no") l2(i)
        write (3,209,advance="no") l3(i)
        write (3,105) ""
    enddo
    
    write (3,206) "P_VECTOR"
    do i=1,panel_num
        write (3,209,advance="no") p1(i)
        write (3,209,advance="no") p2(i)
        write (3,209,advance="no") p3(i)
        write (3,105) ""
    enddo
endif

! writing center points of panels
if (results%cent .eq. 1) then
    write (3,206) "CENTER POINTS"
    do i=1,panel_num
        write (3,209,advance="no") cx(i)
        write (3,209,advance="no") cy(i)
        write (3,209,advance="no") cz(i)
        write (3,105) ""
    enddo
endif

! writing velocities on panels (columns represent cases)
if (results%velo .eq. 1) then
    write (3,206) "VELOCITY"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") v(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing pressure coefficients on panels (columns represent cases)
if (results%pres .eq. 1) then
    write (3,206) "CP"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") cp(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing dipole strengths on panels (columns represent cases)
if (results%doub .eq. 1) then
    write (3,206) "DIPOLE"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") gama(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing source strengths on panels (columns represent cases)
! if requested and if available (doublet/source combination only)
if (results%sorc .eq. 1 .and. method .ne. 1) then
    write (3,206) "SOURCE"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") sigma(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing velocity components on panels (columns represent cases)
if (results%velc .eq. 1) then
    write (3,206) "VELOCITY_COMP"
    do j=1,case_num
        do i=1,panel_num_no_wake
            write (3,209,advance="no") vx(i,j)
            write (3,209,advance="no") vy(i,j)
            write (3,209,advance="no") vz(i,j)
            write (3,105) ""
        enddo
    enddo
endif

! writing static pressure on panels (columns represent cases)
if (results%stat .eq. 1) then
    write (3,206) "P_STATIC"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") p_stat(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing dynamic pressure on panels (columns represent cases)
if (results%dyna .eq. 1) then
    write (3,206) "P_DYNAMIC"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") p_dyna(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing manometer pressure on panels (columns represent cases)
if (results%mano .eq. 1) then
    write (3,206) "P_MANOMETER"
    do i=1,panel_num_no_wake
        do j=1,case_num
            write (3,209,advance="no") p_mano(i,j)
        enddo
        write (3,105) ""
    enddo
endif

! writing stripwise lift coefficients
if (results%cl_strip .eq. 1) then
    write (3,206) "CL_STRIP"
    write (3,206) "X_CENT         Y_CENT         Z_CENT         CL_STRIP"
    do i=1,panel_num-panel_num_no_wake
        do j=1,3
             write (3,209,advance="no") strip_center(i,j)
        end do
        do j=1,case_num
             write (3,209,advance="no") cl_strip(i,j)
        end do
        write (3,105) ""
    enddo
endif

! writing stripwise lift coefficients mutliplied with local chord length
if (results%cl_strip .eq. 1) then
    write (3,206) "CL_STRIP*LOC_CHORD"
    write (3,206) "X_CENT         Y_CENT         Z_CENT         CL_STRIP*LOC_CHORD"
    do i=1,panel_num-panel_num_no_wake
        do j=1,3
             write (3,209,advance="no") strip_center(i,j)
        end do
        do j=1,case_num
             write (3,209,advance="no") cllength_strip(i,j)
        end do
        write (3,105) ""
    enddo
endif
        

write (3,206) "end"
close (3)

! FORMATS ======================================================================
100 format(i4)
101 format(4x,i2)
102 format(6x,i2)
103 format(i2)
104 format(2x,i2)
105 format(a)
106 format(1x,"Date:",1x,i2.2,"/",i2.2,"/",i4,2x,"Time:",1x,i2.2,":",i2.2)
206 format(1x,a)
208 format(1x,i10)
209 format(1x,E14.8)
! ==============================================================================

return
end subroutine subr_results
