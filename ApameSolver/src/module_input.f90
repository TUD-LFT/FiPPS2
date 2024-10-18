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

! file module_input.f90

! This module contains subroutine input that reads apame input file

module module_input

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                                             :: case_num                              ,&
                                                     & method                                ,&
                                                     & collcalc                              ,&
                                                     & velorder                              ,&
                                                     & symmetry                              ,&
                                                     & trim_cl                               ,&
                                                     & wake_mod                              ,&
                                                     & res_req                               ,&
                                                     & vtk_req                               ,&
                                                     & node_num                              ,&
                                                     & panel_num                             ,&
                                                     & panel_num_no_wake                     ,&
                                                     & input_err
real(kind=kind_float)                               :: speed                                 ,&
                                                     & ro                                    ,&
                                                     & p_ref                                 ,&
                                                     & mach                                  ,&
                                                     & wing_span                             ,&
                                                     & mac                                   ,&
                                                     & wing_surf                             ,&
                                                     & origin(3)                             ,&
                                                     & r_stabilizer                          ,&
                                                     & trefftz(3)                            ,&
                                                     & target_cl                             ,&
                                                     & cl_eps                                ,&
                                                     & error                                 ,&
                                                     & coplanarity_angle                     ,&
                                                     & colldist                              ,&
                                                     & farfield
integer, allocatable, dimension(:)                  :: panel_type                            ,&
                                                     & node1,node2,node3,node4               ,&
                                                     & elem1,elem2,elem3,elem4
real(kind=kind_float),    allocatable, dimension(:) :: alfa,beta                             ,&
                                                     & x_node,y_node,z_node
TYPE(results_type)                                  :: results
! ==============================================================================

contains
    
    subroutine input( interactive                                           ,&
                    & file_name                                             ,&
                    & version                                               ,&
                    & subversion                                            )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,            intent(in)                :: interactive
    character(len=*),   intent(in)                :: version,subversion,file_name
    ! PRIVATE ======================================================================
    integer                                       :: i,j                                         ,&
                                                   & IOstatus                                    ,&
                                                   & current_panel_type                          ,&
                                                   & num_lines                                   ,&
                                                   & body_panels_read
    character(len=100)                            :: read_version                                ,&
                                                   & keyword                                     ,&
                                                   & value
    character(len=200)                            :: line
    character(len=200), allocatable, dimension(:) :: file_lines

    ! ==============================================================================

    ! open input file for reading
    call func_message( interactive                                          ,&
                     & "Reading APAME input file...")
    open(unit=1,file=trim(file_name)//".inp",status="old")
    
    ! reading number of lines
    num_lines=0
    do
        read (1,"(a1)",IOSTAT=IOstatus) line
        if (IOstatus .eq. 0) then
            num_lines=num_lines+1
        else
            exit
        endif
    enddo
    rewind(1)

    ! store lines into array of lines
    allocate(file_lines(num_lines))
    do i=1,num_lines
        read (1,"(a)") file_lines(i)
    enddo
    close(1)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! fill default keyword values
    speed=1                                 ! airspeed
    ro=1.225_kind_float                     ! density
    p_ref=101325                            ! pressure
    mach=0                                  ! mach number
    case_num=1                              ! number of cases
    allocate(alfa(1),beta(1))               ! flow angles
    alfa(1)=0
    beta(1)=0
    wing_span=1                             ! wing span
    mac=1                                   ! mean aerodynamic chord
    wing_surf=1                             ! wing surface area
    origin=(/0,0,0/)                        ! origin for moments calculation
    r_stabilizer=0._kind_float              ! distance from origin to center of pressure of the 
                                            ! horizontal stabilizer
    trefftz=(/40,0,0/)                      ! origin of trefftz-plane for wake modification
    method=0                                ! singularity method
    error=0.0000001_kind_float              ! error/small value
    coplanarity_angle=179*pi/180            ! coplanarity angle for curvature correction
    colldist=0.0000001_kind_float           ! collocation point depth
    farfield=5                              ! farfield coefficient
    collcalc=0                              ! type of collocation point calculation
    velorder=1                              ! interpolation method/order for velocity calculations
    symmetry=0                              ! Symmetry to x-z-plane
    trim_cl=0                               ! Trim alfa to target lift coefficient
    target_cl=0._kind_float                 ! target lift coefficient
    cl_eps=1.e-3_kind_float                 ! epsilon used for trimming lift coefficient
    wake_mod=0                              ! Modify wake
    res_req=1                               ! results request
    vtk_req=0                               ! vtk-file request
    results%coef=1                          ! coefficients
    results%forc=0                          ! forces
    results%geom=0                          ! geometry
    results%velo=0                          ! velocities
    results%pres=0                          ! pressure coefficients
    results%cent=0                          ! center points
    results%doub=0                          ! doublet values
    results%sorc=0                          ! source values
    results%velc=0                          ! velocity components
    results%mesh=0                          ! mesh characteristics
    results%stat=0                          ! static pressure
    results%dyna=0                          ! dynamic pressure
    results%mano=0                          ! manometer pressure
    results%cl_strip=0                      ! stripwise lift coefficient
    
    ! assume no error in input module occurred
    input_err=0
    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! reading starting APAME keyword
    if (trim(adjustl(file_lines(1))) .ne. "APAME input file") then
        call func_message( interactive                                      ,&
                         & "ERROR: Not an APAME input file")
        input_err=1
        goto 999
    endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! reading VERSION keyword
    value="0"
    read (file_lines(2),*,IOSTAT=IOstatus) keyword,value
    if (IOstatus .ne. 0) then
        call func_message( interactive                                      ,&
                         & "ERROR: VERSION sintax incorrect")
        input_err=2
        goto 999
    endif

    if (keyword .ne. "VERSION") then
        call func_message( interactive                                      ,&
                         & "ERROR: <VERSION> keyword should appear after <APAME input file> inside input file")
        input_err=1
        goto 999
    endif

    ! reading value for VERSION
    read(value,*) read_version
    if (read_version .ne. trim(version)//"."//trim(subversion)) then
        call func_message( interactive                                      ,&
                         & "ERROR: Wrong input file version detected")
        input_err=1
        goto 999
    endif

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! READING KEYWORDS AND THEIR VALUES
    i=2
    do while (i .lt. num_lines)
        i=i+1
        
        ! get current line and trim for blanks
        line=trim(adjustl(file_lines(i)))
        
        ! check if comment or empty line
        if (line(1:1) .eq. "#" .or. line(1:1) .eq. "") cycle
        
        ! read keyword and value
        read (line,*,IOSTAT=IOstatus) keyword,value
        if (IOstatus .ne. 0) then
            call func_message( interactive                                  ,&
                             & "ERROR: Keyword <"//trim(keyword)//"> does not have a value")
            input_err=2
            goto 999
        endif
        
        ! check for possible keywords
        if (keyword .eq. "AIRSPEED") then
            read(value,*) speed
            if (speed .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: AIRSPEED parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "DENSITY") then
            read(value,*) ro
            if (ro .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: DENSITY parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "PRESSURE") then
            read(value,*) p_ref
            if (p_ref .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: PRESSURE parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "MACH") then
            read(value,*) mach
            if (mach .lt. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: MACH parameter should be equal or greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "CASE_NUM") then
            read(value,*) case_num
            if (case_num .le. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: CASE_NUM parameter should be greater then zero")
                input_err=2
                goto 999
            endif
            
            deallocate(alfa,beta)
            allocate(alfa(case_num),beta(case_num))
            
            ! check that extra 2 lines exist for ALFA and BETA angles
            if (i+2 .gt. num_lines) then
                call func_message( interactive                              ,&
                                 & "ERROR: Unexpected end of file")
                input_err=1
                goto 999
            endif
            
            ! read AOA values
            i=i+1
            line=trim(file_lines(i))
            read (line,*,IOSTAT=IOstatus) alfa
            if (IOstatus .ne. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: CASE_NUM parameter indicates more angles of attack then available")
                input_err=2
                goto 999
            endif
            alfa=alfa*pi/180
            
            ! read BETA values
            i=i+1
            line=trim(file_lines(i))
            read (line,*,IOSTAT=IOstatus) beta
            if (IOstatus .ne. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: CASE_NUM parameter indicates more sideslip angles then available")
                input_err=2
                goto 999
            endif
            beta=beta*pi/180
        elseif (keyword .eq. "WINGSPAN") then
            read(value,*) wing_span
            if (wing_span .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: WINGSPAN parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "MAC") then
            read(value,*) mac
            if (mac .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: MAC parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "SURFACE") then
            read(value,*) wing_surf
            if (wing_surf .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: SURFACE parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "ORIGIN") then
            ! check that extra line exists for origin values
            if (i+1 .gt. num_lines) then
                call func_message( interactive                              ,&
                                 & "ERROR: Unexpected end of file")
                input_err=1
                goto 999
            endif
            i=i+1
            line=trim(file_lines(i))
            read (line,*,IOSTAT=IOstatus) origin
            if (IOstatus .ne. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: Error in declaring (x,y,z) origin coordinates")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "METHOD") then
            read(value,*) method
            if (method .lt. 0 .or. method .gt. 1) then
                call func_message( interactive                              ,&
                                 & "ERROR: METHOD parameter should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "ERROR") then
            read(value,*) error
            if (error .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: ERROR parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "COPL_ANG") then
            read(value,*) coplanarity_angle
            if (error .lt. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: COPL_ANG parameter should be greater or equal to zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "COLLDIST") then
            read(value,*) colldist
            if (colldist .lt. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: COLLDIST parameter should be greater or equal to zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "FARFIELD") then
            read(value,*) farfield
            if (farfield .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: FARFIELD parameter should be greater then zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "COLLCALC") then
            read(value,*) collcalc
            if (collcalc .lt. 0 .or. collcalc .gt. 1) then
                call func_message( interactive                              ,&
                                 & "ERROR: COLLCALC parameter should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "VELORDER") then
            read(value,*) velorder
            if (velorder .lt. 0 .or. velorder .gt. 2) then
                call func_message( interactive                              ,&
                                 & "ERROR: VELORDER parameter should be either 0, 1 or 2")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "SYMMETRY") then
            read(value,*) symmetry
            if (symmetry .lt. 0 .or. symmetry .gt. 1) then
                call func_message( interactive                              ,&
                                 & "ERROR: SYMMETRY parameter should be either 0, or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "TRIM_CL") then
            trim_cl = 1
            ! check that extra line exists for value of target lift coefficient
            if (i+1 .gt. num_lines) then
                call func_message( interactive                              ,&
                                 & "ERROR: Unexpected end of file")
                input_err=1
                goto 999
            endif
            
            ! read desired lift coefficient
            i=i+1
            line=trim(file_lines(i))
            read (line,*,IOSTAT=IOstatus) target_cl
            if (IOstatus .ne. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: Error in declaring target lift coefficient")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "R_STABILIZER") then
            read(value,*,IOSTAT=IOstatus) r_stabilizer
            if (IOstatus .ne. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: Error in declaring distance between horizontal stabilizer and origin")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "CL_EPS") then
            read(value,*) cl_eps
            if (cl_eps .le. 0._kind_float) then
                call func_message( interactive                              ,&
                                 & "ERROR: CL_EPS parameter should be greater zero")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "WAKE_MOD") then
            wake_mod = 1
            ! check that extra line exists for origin values of trefftz-plane
            if (i+1 .gt. num_lines) then
                call func_message( interactive                              ,&
                                 & "ERROR: Unexpected end of file")
                input_err=1
                goto 999
            endif
            i=i+1
            line=trim(file_lines(i))
            read (line,*,IOSTAT=IOstatus) trefftz
            if (IOstatus .ne. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: Error in declaring (x,y,z) origin coordinates of trefftz plane")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RESULTS") then
            read(value,*) res_req
            if (res_req .gt. 1 .or. res_req .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: RESULTS value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "VTK") then
            read(value,*) vtk_req
            if (vtk_req .gt. 1 .or. vtk_req .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_COEF") then
            read(value,*) results%coef
            if (results%coef .gt. 1 .or. results%coef .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_FORC") then
            read(value,*) results%forc
            if (results%forc .gt. 1 .or. results%forc .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_GEOM") then
            read(value,*) results%geom
            if (results%geom .gt. 1 .or. results%geom .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_VELO") then
            read(value,*) results%velo
            if (results%velo .gt. 1 .or. results%velo .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_PRES") then
            read(value,*) results%pres
            if (results%pres .gt. 1 .or. results%pres .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_CENT") then
            read(value,*) results%cent
            if (results%cent .gt. 1 .or. results%cent .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_DOUB") then
            read(value,*) results%doub
            if (results%doub .gt. 1 .or. results%doub .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_SORC") then
            read(value,*) results%sorc
            if (results%sorc .gt. 1 .or. results%sorc .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_VELC") then
            read(value,*) results%velc
            if (results%velc .gt. 1 .or. results%velc .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_MESH") then
            read(value,*) results%mesh
            if (results%mesh .gt. 1 .or. results%mesh .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_STAT") then
            read(value,*) results%stat
            if (results%stat .gt. 1 .or. results%stat .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_DYNA") then
            read(value,*) results%dyna
            if (results%dyna .gt. 1 .or. results%dyna .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_MANO") then
            read(value,*) results%mano
            if (results%mano .gt. 1 .or. results%mano .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "RES_CLST") then
            read(value,*) results%cl_strip
            if (results%cl_strip .gt. 1 .or. results%cl_strip .lt. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: "//trim(keyword)//" value should be either 0 or 1")
                input_err=2
                goto 999
            endif
        elseif (keyword .eq. "NODES") then
            read(value,*) node_num
            if (node_num .le. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: NODES parameter should be greater then 0")
                input_err=2
                goto 999
            endif
            
            ! check that extra lines exist for nodes definition
            if (i+node_num .gt. num_lines) then
                call func_message( interactive                              ,&
                                 & "ERROR: Unexpected end of file... insufficient nodes")
                input_err=1
                goto 999
            endif
            
            ! read nodes
            allocate(x_node(node_num),y_node(node_num),z_node(node_num))
            do j=1,node_num
                i=i+1
                line=trim(file_lines(i))
                read(line,*,IOSTAT=IOstatus) x_node(j),y_node(j),z_node(j)
                if (IOstatus .ne. 0) then
                    call func_message( interactive                          ,&
                                     & "ERROR: Improper node declaration inside input file")
                    input_err=2
                    goto 999
                endif
            enddo
        elseif (keyword .eq. "PANELS") then
            body_panels_read=0
            panel_num_no_wake=0
            read(value,*) panel_num
            if (panel_num .le. 0) then
                call func_message( interactive                              ,&
                                 & "ERROR: PANELS parameter should be greater then 0")
                input_err=2
                goto 999
            endif
            
            ! check that extra lines exist for panels definition
            if (i+panel_num .gt. num_lines) then
                call func_message( interactive                              ,&
                                 & "ERROR: Unexpected end of file... insufficient panels")
                input_err=2
                goto 999
            endif
            
            ! read panels
            allocate( panel_type(panel_num)                                 ,&
                    & node1(panel_num)                                      ,&
                    & node2(panel_num)                                      ,&
                    & node3(panel_num)                                      ,&
                    & node4(panel_num)                                      ,&
                    & elem1(panel_num)                                      ,&
                    & elem2(panel_num)                                      ,&
                    & elem3(panel_num)                                      ,&
                    & elem4(panel_num)                                      )
            do j=1,panel_num
                i=i+1
                line=trim(file_lines(i))
                
                ! read element type
                read (line,*,IOSTAT=IOstatus) current_panel_type
                if (IOstatus .ne. 0) then
                    call func_message( interactive                          ,&
                                     & "ERROR: Improper panel declaration inside input file")
                    input_err=2
                    goto 999
                endif
                
                ! check if body panel defined after already defined wake panel
                if ( current_panel_type .eq. 1                              .or. &
                   & current_panel_type .eq. 2                              .or. &
                   & current_panel_type .eq. 20                             .or. &
                   & current_panel_type .eq. 21                             ) then
                    if (body_panels_read .eq. 1) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Body panels (types 1, 2, 20 and 21) should be defined&
                                         & before wake panels (types 10 and 11)")
                        input_err=2
                        goto 999
                    endif
                endif
                
                if (current_panel_type .eq. 1) then
                    ! quadrilateral body panel
                    
                    ! add to total number of body panels
                    panel_num_no_wake=panel_num_no_wake+1
                    
                    ! read element definition line
                    read (line,*,IOSTAT=IOstatus) panel_type(j)             ,&
                                                & node1(j)                  ,&
                                                & node2(j)                  ,&
                                                & node3(j)                  ,&
                                                & node4(j)                  ,&
                                                & elem1(j)                  ,&
                                                & elem2(j)                  ,&
                                                & elem3(j)                  ,&
                                                & elem4(j)
                    if (IOstatus .ne. 0) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Improper panel declaration inside input file")
                        input_err=2
                        goto 999
                    endif
                    
                    ! check that all node ID's are greater then 0
                    if ( node1(j) .le. 0                                .or. &
                       & node2(j) .le. 0                                .or. &
                       & node3(j) .le. 0                                .or. &
                       & node4(j) .le. 0                                ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Node IDs in panel description should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are greater or equal to 0
                    elseif ( elem1(j) .lt. 0                            .or. &
                           & elem2(j) .lt. 0                            .or. &
                           & elem3(j) .lt. 0                            .or. &
                           & elem4(j) .lt. 0                            ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Element IDs should be greater or equal to 0")
                        input_err=2
                        goto 999
                    ! check that all node ID's are existing
                    elseif ( node1(j) .gt. node_num                     .or. &
                           & node2(j) .gt. node_num                     .or. &
                           & node3(j) .gt. node_num                     .or. &
                           & node4(j) .gt. node_num                     ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced node ID does not exist")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are existing
                    elseif ( elem1(j) .gt. panel_num                    .or. &
                           & elem2(j) .gt. panel_num                    .or. &
                           & elem3(j) .gt. panel_num                    .or. &
                           & elem4(j) .gt. panel_num                    ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced neighboring element ID does not exist")
                        input_err=2
                        goto 999
                    endif
                elseif (current_panel_type .eq. 2) then
                    ! triangular body panel

                    ! add to total number of body panels
                    panel_num_no_wake=panel_num_no_wake+1
                    
                    ! read element definition line
                    read (line,*,IOSTAT=IOstatus) panel_type(j)             ,&
                                                & node1(j)                  ,&
                                                & node2(j)                  ,&
                                                & node3(j)                  ,&
                                                & elem1(j)                  ,&
                                                & elem2(j)                  ,&
                                                & elem3(j)
                    if (IOstatus .ne. 0) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Improper panel declaration inside input file")
                        input_err=2
                        goto 999
                    endif
                    
                    ! check that all node ID's are greater then 0
                    if ( node1(j) .le. 0                                .or. &
                       & node2(j) .le. 0                                .or. &
                       & node3(j) .le. 0                                ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Node IDs in panel description should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are greater or equal to 0
                    elseif ( elem1(j) .lt. 0                            .or. &
                           & elem2(j) .lt. 0                            .or. &
                           & elem3(j) .lt. 0                            ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Element IDs should be greater or equal to 0")
                        input_err=2
                        goto 999
                    ! check that all node ID's are existing
                    elseif ( node1(j) .gt. node_num                     .or. &
                           & node2(j) .gt. node_num                     .or. &
                           & node3(j) .gt. node_num                     ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced node ID does not exist")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are existing
                    elseif ( elem1(j) .gt. panel_num                    .or. &
                           & elem2(j) .gt. panel_num                    .or. &
                           & elem3(j) .gt. panel_num                    ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced neighboring element ID does not exist")
                        input_err=2
                        goto 999
                    endif
                    
                    ! declaring non existing ID's as 0
                    node4(j)=0
                    elem4(j)=0
                elseif (current_panel_type .eq. 20) then
                    ! quadrilateral dummy panel
                    
                    ! add to total number of body panels
                    panel_num_no_wake=panel_num_no_wake+1
                    
                    ! read element definition line
                    read (line,*,IOSTAT=IOstatus) panel_type(j)             ,&
                                                & node1(j)                  ,&
                                                & node2(j)                  ,&
                                                & node3(j)                  ,&
                                                & node4(j)                  ,&
                                                & elem1(j)                  ,&
                                                & elem2(j)                  ,&
                                                & elem3(j)                  ,&
                                                & elem4(j)
                    if (IOstatus .ne. 0) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Improper panel declaration inside input file")
                        input_err=2
                        goto 999
                    endif
                    
                    ! check that all node ID's are greater then 0
                    if ( node1(j) .le. 0                                .or. &
                       & node2(j) .le. 0                                .or. &
                       & node3(j) .le. 0                                .or. &
                       & node4(j) .le. 0                                ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Node IDs in panel description should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are greater or equal to 0
                    elseif ( elem1(j) .lt. 0                            .or. &
                           & elem2(j) .lt. 0                            .or. &
                           & elem3(j) .lt. 0                            .or. &
                           & elem4(j) .lt. 0                            ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Element IDs should be greater or equal to 0")
                        input_err=2
                        goto 999
                    ! check that all node ID's are existing
                    elseif ( node1(j) .gt. node_num                     .or. &
                           & node2(j) .gt. node_num                     .or. &
                           & node3(j) .gt. node_num                     .or. &
                           & node4(j) .gt. node_num                     ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced node ID does not exist")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are existing
                    elseif ( elem1(j) .gt. panel_num                    .or. &
                           & elem2(j) .gt. panel_num                    .or. &
                           & elem3(j) .gt. panel_num                    .or. &
                           & elem4(j) .gt. panel_num                    ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced neighboring element ID does not exist")
                        input_err=2
                        goto 999
                    endif
                elseif (current_panel_type .eq. 21) then
                    ! triangular dummy panel
                    
                    ! add to total number of body panels
                    panel_num_no_wake=panel_num_no_wake+1
                    
                    ! read element definition line
                    read (line,*,IOSTAT=IOstatus) panel_type(j)             ,&
                                                & node1(j)                  ,&
                                                & node2(j)                  ,&
                                                & node3(j)                  ,&
                                                & elem1(j)                  ,&
                                                & elem2(j)                  ,&
                                                & elem3(j)
                    if (IOstatus .ne. 0) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Improper panel declaration inside input file")
                        input_err=2
                        goto 999
                    endif
                    
                    ! check that all node ID's are greater then 0
                    if ( node1(j) .le. 0                                .or. &
                       & node2(j) .le. 0                                .or. &
                       & node3(j) .le. 0                                ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Node IDs in panel description should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are greater or equal to 0
                    elseif ( elem1(j) .lt. 0                            .or. &
                           & elem2(j) .lt. 0                            .or. &
                           & elem3(j) .lt. 0                            ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Element IDs should be greater or equal to 0")
                        input_err=2
                        goto 999
                    ! check that all node ID's are existing
                    elseif ( node1(j) .gt. node_num                     .or. &
                           & node2(j) .gt. node_num                     .or. &

                           & node3(j) .gt. node_num                     ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced node ID does not exist")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are existing
                    elseif ( elem1(j) .gt. panel_num                    .or. &
                           & elem2(j) .gt. panel_num                    .or. &
                           & elem3(j) .gt. panel_num                    ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced neighboring element ID does not exist")
                        input_err=2
                        goto 999
                    endif
                    
                    ! declaring non existing ID's as 0
                    node4(j)=0
                    elem4(j)=0
                elseif (current_panel_type .eq. 10) then
                    ! quadrilateral wake panel
                    
                    ! flag that body panels have been read
                    body_panels_read=1
                    
                    ! read element definition line
                    read (line,*,IOSTAT=IOstatus) panel_type(j)             ,&
                                                & node1(j)                  ,&
                                                & node2(j)                  ,&
                                                & node3(j)                  ,&
                                                & node4(j)                  ,&
                                                & elem1(j)                  ,&
                                                & elem2(j)
                    if (IOstatus .ne. 0) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Improper panel declaration inside input file")
                        input_err=2
                        goto 999
                    endif
                    
                    ! check that all node ID's are greater then 0
                    if ( node1(j) .le. 0                                .or. &
                       & node2(j) .le. 0                                .or. &
                       & node3(j) .le. 0                                .or. &
                       & node4(j) .le. 0                                ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Node IDs in panel description should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that trailing element ID's are greater then 0
                    elseif ( elem1(j) .le. 0                            .or. &
                           & elem2(j) .le. 0                            ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Wake element IDs should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that all node ID's are existing
                    elseif ( node1(j) .gt. node_num                     .or. &
                           & node2(j) .gt. node_num                     .or. &
                           & node3(j) .gt. node_num                     .or. &
                           & node4(j) .gt. node_num                     ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced node ID does not exist")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are existing
                    elseif ( elem1(j) .gt. panel_num                    .or. &
                           & elem2(j) .gt. panel_num                    ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced neighboring element ID does not exist")
                        input_err=2
                        goto 999
                    endif
                    
                    ! declaring non existing ID's as 0
                    elem3(j)=0
                    elem4(j)=0
                elseif (current_panel_type .eq. 11) then
                    ! triangular wake panel
                    
                    ! flag that body panels have been read
                    body_panels_read=1
                    
                    read (line,*,IOSTAT=IOstatus) panel_type(j)             ,&
                                                & node1(j)                  ,&
                                                & node2(j)                  ,&
                                                & node3(j)                  ,&
                                                & elem1(j)                  ,&
                                                & elem2(j)                  ,&
                                                & elem3(j)
                    if (IOstatus .ne. 0) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Improper panel declaration inside input file")
                        input_err=2
                        goto 999
                    endif
                    
                    ! check that all node ID's are greater then 0
                    if ( node1(j) .le. 0                                .or. &
                       & node2(j) .le. 0                                .or. &
                       & node3(j) .le. 0                                ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Node IDs in panel description should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that trailing element ID's are greater then 0
                    elseif ( elem1(j) .le. 0                            .or. &
                           & elem2(j) .le. 0                            ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Wake element IDs should be greater then 0")
                        input_err=2
                        goto 999
                    ! check that all node ID's are existing
                    elseif ( node1(j) .gt. node_num                     .or. &
                           & node2(j) .gt. node_num                     .or. &
                           & node3(j) .gt. node_num                     ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced node ID does not exist")
                        input_err=2
                        goto 999
                    ! check that all neighboring element ID's are existing
                    elseif ( elem1(j) .gt. panel_num                    .or. &
                           & elem2(j) .gt. panel_num                    ) then
                        call func_message( interactive                      ,&
                                         & "ERROR: Referenced neighboring element ID does not exist")
                        input_err=2
                        goto 999
                    endif
                    
                    ! declaring non existing ID's as 0
                    node4(j)=0
                    elem3(j)=0
                    elem4(j)=0
                else
                    call func_message( interactive                          ,&
                                     & "ERROR: Panel type should be either 1, 2, 10, 11, 20 or 21")
                    input_err=2
                    goto 999
                endif
            enddo
        else
            ! unknown keyword found
            call func_message( interactive                                  ,&
                             & "ERROR: Unknown keyword <"//trim(keyword)//"> found")
            input_err=2
            goto 999
        endif
    enddo
    
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

999 continue

    if (input_err .eq. 2) then
        call func_new_line(interactive)
        write (line,'(I10.10)') i
        call func_message( interactive                                      ,&
                         & "Error at line "//trim(line))
    endif
    
    return
    end subroutine input

end module module_input

