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

! file module_velocity.f90

! This module calculates induced velocities from doublet strengths
! Inverse Distance Weighting method used in nodal interpolation method
! developed by Julien Chaussee

module module_velocity

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                                            :: velo_err
real(kind=kind_float), allocatable, dimension(:,:) :: v,vx,vy,vz
! ==============================================================================

contains

! SUBROUTINE VELOCITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine velocities( gama                                             ,&
                         & node_num                                         ,&
                         & panel_num                                        ,&
                         & panel_num_no_wake                                ,&
                         & case_num                                         ,&
                         & velorder                                         ,&
                         & method                                           ,&
                         & symmetry                                         ,&
                         & coplanarity_angle                                ,&
                         & panel_type                                       ,&
                         & interactive                                      ,&
                         & cx,cy,cz                                         ,&
                         & x,y,z                                            ,&
                         & l,p,n                                            ,&
                         & speed                                            ,&
                         & speed_x,speed_y,speed_z                          ,&
                         & node1,node2,node3,node4                          ,&
                         & elem1,elem2,elem3,elem4                          ,&
                         & error                                            )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                                intent(in)                  :: node_num          ,&
                                                                         & panel_num         ,&
                                                                         & panel_num_no_wake ,&
                                                                         & interactive       ,&
                                                                         & case_num          ,&
                                                                         & velorder          ,&
                                                                         & method            ,&
                                                                         & symmetry
    real(kind=kind_float),                                   intent(in) :: speed             ,&
                                                                         & coplanarity_angle ,&
                                                                         & error
    integer, dimension(panel_num),          intent(in)                  :: panel_type        ,&
                                                                         & node1,node2       ,&
                                                                         & node3,node4
    integer, dimension(panel_num),          intent(inout)               :: elem1,elem2       ,&
                                                                         & elem3,elem4
    real(kind=kind_float),    dimension(panel_num),          intent(in) :: cx,cy,cz
    real(kind=kind_float),    dimension(node_num),           intent(in) :: x,y,z
    real(kind=kind_float),    dimension(case_num),           intent(in) :: speed_x,speed_y,speed_z
    real(kind=kind_float),    dimension(panel_num,3),        intent(in) :: l,p,n
    real(kind=kind_float),    dimension(panel_num_no_wake,case_num), intent(in) :: gama
    ! PRIVATE ======================================================================
    integer                                                             :: i,j               ,&
                                                                         & v_err             ,&
                                                                         & interp_error      ,&
                                                                         & alloc_stat
    integer, dimension(panel_num)                                       :: elem_sym, temp
    integer, allocatable, dimension(:)                                  :: trailing_panel
    real(kind=kind_float)                                               :: gl,gp
    real(kind=kind_float),    allocatable, dimension(:)                 :: ql,qp
    real(kind=kind_float),    allocatable, dimension(:,:)               :: gama_nodal        ,&
                                                                         & node_weight_sum   ,&
                                                                         & weights
    ! ==============================================================================

    ! presume no interpolation errors occured
    interp_error=0
    
    ! presume no velocities subroutine allocation errors
    velo_err=0

    allocate (  v(panel_num_no_wake,case_num)                               ,&
             & vx(panel_num_no_wake,case_num)                               ,&
             & vy(panel_num_no_wake,case_num)                               ,&
             & vz(panel_num_no_wake,case_num)                               ,&
             & stat=alloc_stat                                              )
    if (alloc_stat .ne. 0) then
        call func_message( interactive                                      ,&
                         & "    ERROR: Not enough memory for allocating velocity data fields")
        velo_err=1
        goto 999
    endif

    elem_sym(:) = 0
    ! if symmetry about x-z-plane exists,
    ! find neighboring panels because of symmetry
    if (symmetry .eq. 1) then
        do i=1,panel_num_no_wake
            ! check, if panel is at symmetry plane
            if ((abs(y(node1(i))) .lt. error) .or. (abs(y(node3(i))) .lt. error)) then
                ! if three neighbors exist
                if (elem4(i) .eq. 0 .and. elem3(i) .ne. 0) then
                    ! check, which neighboring panel is not at symmetry plane
                    if (.not. ((abs(y(node1(elem1(i)))) .lt. error) .or. (abs(y(node3(elem1(i)))) .lt. error))) then
                        elem_sym(i) = elem1(i)
                    else if (.not. ((abs(y(node1(elem2(i)))) .lt. error) .or. (abs(y(node3(elem2(i)))) .lt. error))) then
                        elem_sym(i) = elem2(i)
                    else if (.not. ((abs(y(node1(elem3(i)))) .lt. error) .or. (abs(y(node3(elem3(i)))) .lt. error))) then
                        elem_sym(i) = elem3(i)
                    end if
                ! if two neighbors exist
                else if (elem3(i) .eq. 0 .and. elem2(i) .ne. 0) then
                    ! check, which neighboring panel is not at symmetry plane
                    if (.not. ((abs(y(node1(elem1(i)))) .lt. error) .or. (abs(y(node3(elem1(i)))) .lt. error))) then
                        elem_sym(i) = elem1(i)
                    else if (.not. ((abs(y(node1(elem2(i)))) .lt. error) .or. (abs(y(node3(elem2(i)))) .lt. error))) then
                        elem_sym(i) = elem2(i)
                    end if
                ! if one neighbor exists
                else if (elem2(i) .eq. 0 .and. elem1(i) .ne. 0) then
                    ! check, which neighboring panel is not at symmetry plane
                    if (.not. ((abs(y(node1(elem1(i)))) .lt. error) .or. (abs(y(node3(elem1(i)))) .lt. error))) then
                        elem_sym(i) = elem1(i)
                    end if
                end if
            end if
        end do
    end if

    temp(:) = 0
    do i=1,panel_num_no_wake
        if (elem_sym(i) .ne. 0) temp(i) = 1
    end do    

    if (velorder .eq. 0) then
        ! nodal interpolation method
        
        ! allocating variables
        allocate ( weights(panel_num,4)                                     ,&
                 & trailing_panel(panel_num_no_wake)                        ,&
                 & gama_nodal(node_num,case_num)                            ,&
                 & node_weight_sum(node_num,case_num)                       ,&
                 & ql(case_num)                                             ,&
                 & qp(case_num)                                             ,&
                 & stat=alloc_stat                                          )
        if (alloc_stat .ne. 0) then
            call func_message( interactive                                  ,&
                             & "    ERROR: Not enough memory for allocating nodal interpolation variables")
            velo_err=1
            goto 999
        endif
        
        ! calculating weights
        call weight_factors( weights                                        ,&
                           & node1                                          ,&
                           & node2                                          ,&
                           & node3                                          ,&
                           & node4                                          ,&
                           & x                                              ,&
                           & y                                              ,&
                           & z                                              ,&
                           & cx                                             ,&
                           & cy                                             ,&
                           & cz                                             ,&
                           & node_num                                       ,&
                           & panel_num                                      ,&
                           & panel_type)
        
        ! reset values to zero
        gama_nodal(:,:)=0._kind_float
        node_weight_sum(:,:)=0._kind_float
        trailing_panel(:)=0
        
        ! for each panel
        do i=1,panel_num_no_wake
            do j=1,case_num
                gama_nodal(node1(i),j)=gama_nodal(node1(i),j) + gama(i,j)*weights(i,1)
                gama_nodal(node2(i),j)=gama_nodal(node2(i),j) + gama(i,j)*weights(i,2)
                gama_nodal(node3(i),j)=gama_nodal(node3(i),j) + gama(i,j)*weights(i,3)
                node_weight_sum(node1(i),j)=node_weight_sum(node1(i),j) + weights(i,1)
                node_weight_sum(node2(i),j)=node_weight_sum(node2(i),j) + weights(i,2)
                node_weight_sum(node3(i),j)=node_weight_sum(node3(i),j) + weights(i,3)
                if ( panel_type(i) .eq. 1   .or.&
                   & panel_type(i) .eq. 20  )then
                    gama_nodal(node4(i),j)=gama_nodal(node4(i),j) + gama(i,j)*weights(i,4)
                    node_weight_sum(node4(i),j)=node_weight_sum(node4(i),j) + weights(i,4)
                endif
            enddo
        enddo

        ! divide summed influences with weight sum
        do i=1,node_num
            do j=1,case_num
                if (node_weight_sum(i,j) .ne. 0) then
                    gama_nodal(i,j)=gama_nodal(i,j)/node_weight_sum(i,j);
                else
                    gama_nodal(i,j)=0;
                endif
            enddo
        enddo
        
        ! flag trailing panels
        do i=panel_num_no_wake+1,panel_num
            trailing_panel(elem1(i))=1
            trailing_panel(elem2(i))=1
        enddo
        
        ! now interpolate linear surfaces and find velocities
        do i=1,panel_num_no_wake
            if (trailing_panel(i) .eq. 1) then
                ! for trailing panels assume error occured so that free stream velocity is obtained
                v_err=1
            else
                if (panel_type(i) .eq. 1) then
                    do j=1,case_num
                        call linear_4( ql(j),qp(j)                                                      ,&
                                     & (/x(node1(i)),x(node2(i)),x(node3(i)),x(node4(i))/)              ,&
                                     & (/y(node1(i)),y(node2(i)),y(node3(i)),y(node4(i))/)              ,&
                                     & (/z(node1(i)),z(node2(i)),z(node3(i)),z(node4(i))/)              ,&
                                     & l(i,:),p(i,:)                                                    ,&
                                     & (/n(i,:),n(i,:),n(i,:),n(i,:)/)                                  ,&
                                     & v_err                                                            ,&
                                     & 0.0_kind_float                                                              ,&
                                     & (/ gama_nodal(node1(i),j)                                        ,&
                                     &    gama_nodal(node2(i),j)                                        ,&
                                     &    gama_nodal(node3(i),j)                                        ,&
                                     &    gama_nodal(node4(i),j) /)                                     )
                    enddo
                elseif (panel_type(i) .eq. 2) then
                    do j=1,case_num
                        call linear_3( ql(j),qp(j)                                                      ,&
                                     & (/x(node1(i)),x(node2(i)),x(node3(i))/)                          ,&
                                     & (/y(node1(i)),y(node2(i)),y(node3(i))/)                          ,&
                                     & (/z(node1(i)),z(node2(i)),z(node3(i))/)                          ,&
                                     & l(i,:),p(i,:)                                                    ,&
                                     & (/n(i,:),n(i,:),n(i,:)/)                                         ,&
                                     & v_err                                                            ,&
                                     & 0.0_kind_float                                                              ,&
                                     & (/ gama_nodal(node1(i),j)                                        ,&
                                     &    gama_nodal(node2(i),j)                                        ,&
                                     &    gama_nodal(node3(i),j) /)                                     )
                    enddo
                endif
            endif
            
            ! check for singularities
            if (v_err .eq. 1) then
                interp_error=1
                ! define free stream speed magnitude
                if (method .eq. 0) then
                    ! in case source/doublet combination is used, perturbations are zero
                    ql(:)=0._kind_float
                    qp(:)=0._kind_float
                elseif (method .eq. 1) then
                    ! for doublet only, this combination will be closest to free stream velocity
                    ql(:)=sqrt(speed)
                    qp(:)=ql(:)
                endif
            endif
            
            do j=1,case_num
                ! computing v,vx,vy,vz from local induced velocities ql,qp
                if (method .eq. 0) then         ! source/doublet combination (velocity perturbations)
                    ! calculating free stream components
                    gl=sum(l(i,:)*(/speed_x(j),speed_y(j),speed_z(j)/))
                    gp=sum(p(i,:)*(/speed_x(j),speed_y(j),speed_z(j)/))
                    
                    ! calculating velocity components (substituting from free stream components)
                    vx(i,j)=(gl - ql(j))*l(i,1) + (gp - qp(j))*p(i,1)
                    vy(i,j)=(gl - ql(j))*l(i,2) + (gp - qp(j))*p(i,2)
                    vz(i,j)=(gl - ql(j))*l(i,3) + (gp - qp(j))*p(i,3)
                    
                    ! calculating total velocity
                    v(i,j)=sqrt(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
                elseif (method .eq. 1) then     ! doublets only (total velocity)
                    ! calculating velocity components
                    vx(i,j)=-ql(j)*l(i,1) - qp(j)*p(i,1)
                    vy(i,j)=-ql(j)*l(i,2) - qp(j)*p(i,2)
                    vz(i,j)=-ql(j)*l(i,3) - qp(j)*p(i,3)
                    
                    ! calculating total velocity
                    v(i,j)=sqrt(ql(j)**2 + qp(j)**2)
                endif
            enddo
        enddo
    else
        ! nodal interpolation is not requested so continue surface interpolation method
        allocate ( ql(case_num)                                             ,&
                 & qp(case_num)                                             ,&
                 & stat=alloc_stat                                          )
        if (alloc_stat .ne. 0) then
            call func_message( interactive                                          ,&
                 & "    ERROR: Not enough memory for allocating local velocities")
            velo_err=1
            goto 999
        endif
        
        do i=1,panel_num_no_wake
            ! avoid dummy panels
            if (panel_type(i) .ne. 20 .and. panel_type(i) .ne. 21) then
                ! if four neighbours exist
                if (elem4(i) .ne. 0) then
                    ! if second order interpolation is requested
                    if (velorder .eq. 2) then
                        do j=1,case_num
                        call quadratic_5( ql(j),qp(j)                                                       ,&
                                        & (/cx(i),cx(elem1(i)),cx(elem2(i)),cx(elem3(i)),cx(elem4(i))/)     ,&
                                        & (/cy(i),cy(elem1(i)),cy(elem2(i)),cy(elem3(i)),cy(elem4(i))/)     ,&
                                        & (/cz(i),cz(elem1(i)),cz(elem2(i)),cz(elem3(i)),cz(elem4(i))/)     ,&
                                        & l(i,:),p(i,:)                                                     ,&
                                        & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem3(i),:),n(elem4(i),:)/),&
                                        & v_err                                                             ,&
                                        & coplanarity_angle                                                 ,&
                                        & (/gama(i,j)                                                       ,&
                                        &   gama(elem1(i),j)                                                ,&
                                        &   gama(elem2(i),j)                                                ,&
                                        &   gama(elem3(i),j)                                                ,&
                                        &   gama(elem4(i),j)/)                                              )
                        enddo
                    
                    ! if linear interpolation is requested
                    elseif (velorder .eq. 1) then
                        do j=1,case_num
                        call linear_5( ql(j),qp(j)                                                          ,&
                                     & (/cx(i),cx(elem1(i)),cx(elem2(i)),cx(elem3(i)),cx(elem4(i))/)        ,&
                                     & (/cy(i),cy(elem1(i)),cy(elem2(i)),cy(elem3(i)),cy(elem4(i))/)        ,&
                                     & (/cz(i),cz(elem1(i)),cz(elem2(i)),cz(elem3(i)),cz(elem4(i))/)        ,&
                                     & l(i,:),p(i,:)                                                        ,&
                                     & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem3(i),:),n(elem4(i),:)/)   ,&
                                     & v_err                                                                ,&
                                     & coplanarity_angle                                                    ,&
                                     & (/gama(i,j)                                                          ,&
                                     &   gama(elem1(i),j)                                                   ,&
                                     &   gama(elem2(i),j)                                                   ,&
                                     &   gama(elem3(i),j)                                                   ,&
                                     &   gama(elem4(i),j)/)                                                 )
                        enddo
                    endif
                
                ! if three neighbours exist
                elseif (elem3(i) .ne. 0) then
                    ! if no symmetry must be considered 
                    if (elem_sym(i) .eq. 0) then
                        ! if second order interpolation is requested
                        if (velorder .eq. 2) then
                            do j=1,case_num
                            call quadratic_4( ql(j),qp(j)                                           ,&
                                            & (/cx(i),cx(elem1(i)),cx(elem2(i)),cx(elem3(i))/)      ,&
                                            & (/cy(i),cy(elem1(i)),cy(elem2(i)),cy(elem3(i))/)      ,&
                                            & (/cz(i),cz(elem1(i)),cz(elem2(i)),cz(elem3(i))/)      ,&
                                            & l(i,:),p(i,:)                                         ,&
                                            & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem3(i),:)/)  ,&
                                            & v_err                                                 ,&
                                            & coplanarity_angle                                     ,&
                                            & (/gama(i,j)                                           ,&
                                            &   gama(elem1(i),j)                                    ,&
                                            &   gama(elem2(i),j)                                    ,&
                                            &   gama(elem3(i),j)/)                                  )
                            enddo
                        
                        ! if linear interpolation is requested
                        elseif (velorder .eq. 1) then
                            do j=1,case_num
                            call linear_4( ql(j),qp(j)                                          ,&
                                        & (/cx(i),cx(elem1(i)),cx(elem2(i)),cx(elem3(i))/)     ,&
                                        & (/cy(i),cy(elem1(i)),cy(elem2(i)),cy(elem3(i))/)     ,&
                                        & (/cz(i),cz(elem1(i)),cz(elem2(i)),cz(elem3(i))/)     ,&
                                        & l(i,:),p(i,:)                                        ,&
                                        & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem3(i),:)/) ,&
                                        & v_err                                                ,&
                                        & coplanarity_angle                                    ,&
                                        & (/gama(i,j)                                          ,&
                                        &   gama(elem1(i),j)                                   ,&
                                        &   gama(elem2(i),j)                                   ,&
                                        &   gama(elem3(i),j)/)                                 )
                            enddo
                        endif
                    ! if symmetry must be considered
                    else
                        ! if second order interpolation is requested
                        if (velorder .eq. 2) then
                            do j=1,case_num
                            call quadratic_5( ql(j),qp(j)                                                      ,&
                                        & (/cx(i),cx(elem1(i)),cx(elem2(i)),cx(elem3(i)), cx(elem_sym(i))/)    ,&
                                        & (/cy(i),cy(elem1(i)),cy(elem2(i)),cy(elem3(i)),-cy(elem_sym(i))/)    ,&
                                        & (/cz(i),cz(elem1(i)),cz(elem2(i)),cz(elem3(i)), cz(elem_sym(i))/)    ,&
                                        & l(i,:),p(i,:)                                                        ,&
                                        & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem3(i),:),n(elem_sym(i),:)/),&
                                        & v_err                                                                ,&
                                        & coplanarity_angle                                                    ,&
                                        & (/gama(i,j)                                                          ,&
                                        &   gama(elem1(i),j)                                                   ,&
                                        &   gama(elem2(i),j)                                                   ,&
                                        &   gama(elem3(i),j)                                                   ,&
                                        &   gama(elem_sym(i),j)/)                                              )
                            enddo
                        
                        ! if linear interpolation is requested
                        elseif (velorder .eq. 1) then
                            do j=1,case_num
                            call linear_5( ql(j),qp(j)                                                          ,&
                                     & (/cx(i),cx(elem1(i)),cx(elem2(i)),cx(elem3(i)), cx(elem_sym(i))/)        ,&
                                     & (/cy(i),cy(elem1(i)),cy(elem2(i)),cy(elem3(i)),-cy(elem_sym(i))/)        ,&
                                     & (/cz(i),cz(elem1(i)),cz(elem2(i)),cz(elem3(i)), cz(elem_sym(i))/)        ,&
                                     & l(i,:),p(i,:)                                                            ,&
                                     & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem3(i),:),n(elem_sym(i),:)/)    ,&
                                     & v_err                                                                    ,&
                                     & coplanarity_angle                                                        ,&
                                     & (/gama(i,j)                                                              ,&
                                     &   gama(elem1(i),j)                                                       ,&
                                     &   gama(elem2(i),j)                                                       ,&
                                     &   gama(elem3(i),j)                                                       ,&
                                     &   gama(elem_sym(i),j)/)                                                  )
                            enddo
                        endif
                    endif
                
                ! if two neighbours exist
                elseif (elem2(i) .ne. 0) then
                    ! if no symmetry must be considered 
                    if (elem_sym(i) .eq. 0) then
                        ! only linear interpolation possible
                        do j=1,case_num
                        call linear_3( ql(j),qp(j)                                 ,&
                                    & (/cx(i),cx(elem1(i)),cx(elem2(i))/)          ,&
                                    & (/cy(i),cy(elem1(i)),cy(elem2(i))/)          ,&
                                    & (/cz(i),cz(elem1(i)),cz(elem2(i))/)          ,&
                                    & l(i,:),p(i,:)                                ,&
                                    & (/n(i,:),n(elem1(i),:),n(elem2(i),:)/)       ,&
                                    & v_err                                        ,&
                                    & coplanarity_angle                            ,&
                                    & (/gama(i,j)                                  ,&
                                    &   gama(elem1(i),j)                           ,&
                                    &   gama(elem2(i),j)/)                         )
                        enddo
                    ! if symmetry must be considered
                    else
                        ! if second order interpolation is requested
                        if (velorder .eq. 2) then
                            do j=1,case_num
                            call quadratic_4( ql(j),qp(j)                                            ,&
                                            & (/cx(i),cx(elem1(i)),cx(elem2(i)), cx(elem_sym(i))/)   ,&
                                            & (/cy(i),cy(elem1(i)),cy(elem2(i)),-cy(elem_sym(i))/)   ,&
                                            & (/cz(i),cz(elem1(i)),cz(elem2(i)), cz(elem_sym(i))/)   ,&
                                            & l(i,:),p(i,:)                                          ,&
                                            & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem_sym(i),:)/),&
                                            & v_err                                                  ,&
                                            & coplanarity_angle                                      ,&
                                            & (/gama(i,j)                                            ,&
                                            &   gama(elem1(i),j)                                     ,&
                                            &   gama(elem2(i),j)                                     ,&
                                            &   gama(elem_sym(i),j)/)                                )
                            enddo
                        
                        ! if linear interpolation is requested
                        elseif (velorder .eq. 1) then
                            do j=1,case_num
                            call linear_4( ql(j),qp(j)                                           ,&
                                        & (/cx(i),cx(elem1(i)),cx(elem2(i)), cx(elem_sym(i))/)   ,&
                                        & (/cy(i),cy(elem1(i)),cy(elem2(i)),-cy(elem_sym(i))/)   ,&
                                        & (/cz(i),cz(elem1(i)),cz(elem2(i)), cz(elem_sym(i))/)   ,&
                                        & l(i,:),p(i,:)                                          ,&
                                        & (/n(i,:),n(elem1(i),:),n(elem2(i),:),n(elem_sym(i),:)/),&
                                        & v_err                                                  ,&
                                        & coplanarity_angle                                      ,&
                                        & (/gama(i,j)                                            ,&
                                        &   gama(elem1(i),j)                                     ,&
                                        &   gama(elem2(i),j)                                     ,&
                                        &   gama(elem_sym(i),j)/)                                )
                            enddo
                        endif
                    endif
                
                ! if one or no neighbors exist
                else
                    ! if no symmetry must be considered 
                    if (elem_sym(i) .eq. 0) then
                        ! declare freestream velocity
                        v_err=0
                        if (method .eq. 0) then
                            ! in case source/doublet combination is used, perturbations are zero
                            ql(:)=0._kind_float
                            qp(:)=0._kind_float
                        elseif (method .eq. 1) then
                            ! for doublet only, this combination will be closest to free stream velocity
                            ql(:)=sqrt(speed)
                            qp(:)=ql(:)
                        endif
                    ! if symmetry must be considered
                    else
                        ! only linear interpolation possible
                        do j=1,case_num
                        call linear_3( ql(j),qp(j)                                 ,&
                                    & (/cx(i),cx(elem1(i)), cx(elem_sym(i))/)      ,&
                                    & (/cy(i),cy(elem1(i)),-cy(elem_sym(i))/)      ,&
                                    & (/cz(i),cz(elem1(i)), cz(elem_sym(i))/)      ,&
                                    & l(i,:),p(i,:)                                ,&
                                    & (/n(i,:),n(elem1(i),:),n(elem_sym(i),:)/)    ,&
                                    & v_err                                        ,&
                                    & coplanarity_angle                            ,&
                                    & (/gama(i,j)                                  ,&
                                    &   gama(elem1(i),j)                           ,&
                                    &   gama(elem_sym(i),j)/)                      )
                        enddo
                    endif
                endif
                
                ! check for issues in interpolation
                if (v_err .eq. 1) then
                    interp_error=1
                    ! singularity detected so define free stream speed magnitude
                    if (method .eq. 0) then
                        ! in case source/doublet combination is used, perturbations are zero
                        ql(:)=0._kind_float
                        qp(:)=0._kind_float
                    elseif (method .eq. 1) then
                        ! for doublet only, this combination will be closest to free stream velocity
                        ql(:)=sqrt(speed)
                        qp(:)=ql(:)
                    endif
                endif
                
                do j=1,case_num
                    ! computing v,vx,vy,vz from local induced velocities ql,qp
                    if (method .eq. 0) then         ! source/doublet combination (velocity perturbations)
                        ! calculating free stream components
                        gl=sum(l(i,:)*(/speed_x(j),speed_y(j),speed_z(j)/))
                        gp=sum(p(i,:)*(/speed_x(j),speed_y(j),speed_z(j)/))
                        
                        ! calculating velocity components (substituting from free stream components)
                        vx(i,j)=(gl - ql(j))*l(i,1) + (gp - qp(j))*p(i,1)
                        vy(i,j)=(gl - ql(j))*l(i,2) + (gp - qp(j))*p(i,2)
                        vz(i,j)=(gl - ql(j))*l(i,3) + (gp - qp(j))*p(i,3)
                        
                        ! calculating total velocity
                        v(i,j)=sqrt(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)
                    elseif (method .eq. 1) then     ! doublets only (total velocity)
                        ! calculating velocity components
                        vx(i,j)=- ql(j)*l(i,1) - qp(j)*p(i,1)
                        vy(i,j)=- ql(j)*l(i,2) - qp(j)*p(i,2)
                        vz(i,j)=- ql(j)*l(i,3) - qp(j)*p(i,3)
                        
                        ! calculating total velocity
                        v(i,j)=sqrt(ql(j)**2 + qp(j)**2)
                    endif
                enddo
            endif
        enddo

        if (interp_error .eq. 1) then
            call func_message( interactive                                          ,&
                & "    WARNING: Velocity interpolation singularities detected, check mesh... ")
        endif
    endif
    
999 continue
    
    return
    end subroutine velocities
    
! SUBROUTINE QUADRATIC_5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! This function calculates induced velocities (perturbations or total,
    ! depending on method) based on given five dipole values (gama), their
    ! locations (x,y,z) and local unit vectors (l and p).
    ! Method is based on interpolating 5 points into quadratic surface
    ! with equation gama=A*x + B*y + C*x^2 + D*y^2 + E.
    ! Induced velocities in l and p directions will be coefficients A and B.
    
    subroutine quadratic_5( ql,qp                                           ,&
                          & x,y,z                                           ,&
                          & l,p                                             ,&
                          & n                                               ,&
                          & v_err                                           ,&
                          & coplanarity_angle                               ,&
                          & gama                                            )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                       intent(in)  :: coplanarity_angle
    real(kind=kind_float),    dimension(5),      intent(in)  :: x,y,z,gama
    real(kind=kind_float),    dimension(3),      intent(in)  :: l,p
    real(kind=kind_float),    dimension(3,5),    intent(in)  :: n
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                       intent(out) :: ql,qp
    ! INTENTS INOUT ================================================================
    integer,                    intent(out) :: v_err
    ! PRIVATE ======================================================================
    integer                                                  :: i,info
    integer, dimension(5,5)                                  :: ipiv
    real(kind=kind_float)                                    :: correcting_factor
    real(kind=kind_float),    dimension(3)                   :: radius_vector
    real(kind=kind_float),    dimension(5,5)                 :: T
    ! ==============================================================================

    ! finding point coordinates in local coordinate system and corresponding
    ! "T" coefficients.
    T(1,:)=(/0._kind_float, 0._kind_float, 0._kind_float, 0._kind_float, 1._kind_float/)
    do i=2,5    ! for each of five gama points
        ! find radius vector of "i" point relative to first "1" which is origin of
        ! local coordinate system
        radius_vector=(/x(i)-x(1), y(i)-y(1), z(i)-z(1)/)
        
        ! calculate correction factor that will take into account geometry curvature
        call curvature_correction( correcting_factor                        ,&
                                 & coplanarity_angle                        ,&
                                 & radius_vector                            ,&
                                 & n(:,1)                                   ,&
                                 & n(:,i)                                   )
        
        ! transforming this vector to local coordinate system
        T(i,1)=sum(radius_vector*l)*correcting_factor   ! "l" coordinate in local coordinate system (coeff A)
        T(i,2)=sum(radius_vector*p)*correcting_factor   ! "p" coordinate in local coordinate system (coeff B)
        T(i,3)=T(i,1)**2                              ! "l" squared (coeff C)
        T(i,4)=T(i,2)**2                              ! "p" squared (coeff D)
        T(i,5)=1._kind_float                          ! (coeff E)
    enddo

    ! solving system of equations and returning coefficients in gama vector
    if (kind_float .eq. kind_float_l) then
        ! single precision
        call sgesv( 5                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 5                                                           ,&
                  & ipiv                                                        ,&
                  & gama                                                        ,&
                  & 5                                                           ,&
                  & info                                                        )
    else
        ! double precision
        call dgesv( 5                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 5                                                           ,&
                  & ipiv                                                        ,&
                  & gama                                                        ,&
                  & 5                                                           ,&
                  & info                                                        )
    end if
    
    
    if (info .eq. 0) then
        ! localy induced velocities ql,qp
        ql=gama(1)    ! A coefficient
        qp=gama(2)    ! B coefficient
        v_err=0
    else
        ! problem with interpolation occured
        ql=0._kind_float
        qp=0._kind_float
        v_err=1
    endif

    return
    end subroutine quadratic_5
    
! SUBROUTINE QUADRATIC_4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! This function calculates induced velocities (perturbations or total,
    ! depending on method) based on given four dipole values (gama), their
    ! locations (x,y,z) and local unit vectors (l and p).
    ! Method is based on interpolating 4 points into quadratic surface
    ! with equation gama=A*x + B*y + C*(x^2+y^2) + D.
    ! Induced velocities in l and p directions will be coefficients A and B.
    
    subroutine quadratic_4( ql,qp                                           ,&
                          & x,y,z                                           ,&
                          & l,p                                             ,&
                          & n                                               ,&
                          & v_err                                           ,&
                          & coplanarity_angle                               ,&
                          & gama                                            )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                       intent(in)  :: coplanarity_angle
    real(kind=kind_float),    dimension(4),      intent(in)  :: x,y,z,gama
    real(kind=kind_float),    dimension(3),      intent(in)  :: l,p
    real(kind=kind_float),    dimension(3,4),    intent(in)  :: n
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                       intent(out) :: ql,qp
    ! INTENTS INOUT ================================================================
    integer,                    intent(out)                  :: v_err
    ! PRIVATE ======================================================================
    integer                                 :: i,info
    integer, dimension(4,4)                 :: ipiv
    real(kind=kind_float)                                    :: correcting_factor
    real(kind=kind_float),    dimension(3)                   :: radius_vector
    real(kind=kind_float),    dimension(4,4)                 :: T
    ! ==============================================================================

    ! finding point coordinates in local coordinate system and corresponding
    ! "T" coefficients.
    T(1,:)=(/0._kind_float, 0._kind_float, 0._kind_float, 1._kind_float/)
    do i=2,4    ! for each of four gama points
        ! find radius vector of "i" point relative to first "1" which is origin of
        ! local coordinate system
        radius_vector=(/x(i)-x(1), y(i)-y(1), z(i)-z(1)/)
        
        ! calculate correction factor that will take into account geometry curvature
        call curvature_correction( correcting_factor                        ,&
                                 & coplanarity_angle                        ,&
                                 & radius_vector                            ,&
                                 & n(:,1)                                   ,&
                                 & n(:,i)                                   )
        
        ! transforming this vector to local coordinate system
        T(i,1)=sum(radius_vector*l)       ! "l" coordinate in local coordinate system (coeff A)
        T(i,2)=sum(radius_vector*p)       ! "p" coordinate in local coordinate system (coeff B)
        T(i,3)=T(i,1)**2 + T(i,2)**2      ! "l**2" + "p**2" (coeff C)
        T(i,4)=1._kind_float              ! (coeff D)
    enddo
    
    ! solving system of equations and returning coefficients in gama vector
    if (kind_float .eq. kind_float_l) then
        ! single precision
        call sgesv( 4                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 4                                                           ,&
                  & ipiv                                                        ,&
                  & gama                                                        ,&
                  & 4                                                           ,&
                  & info                                                        )
    else
        ! double precision
        call dgesv( 4                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 4                                                           ,&
                  & ipiv                                                        ,&
                  & gama                                                        ,&
                  & 4                                                           ,&
                  & info                                                        )
    end if

    if (info .eq. 0) then
        ! localy induced velocities ql,qp
        ql=gama(1)    ! A coefficient
        qp=gama(2)    ! B coefficient
        v_err=0
    else
        ! problem with interpolation occured
        ql=0._kind_float
        qp=0._kind_float
        v_err=1
    endif

    return
    end subroutine quadratic_4
    
! SUBROUTINE LINEAR_5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! This function calculates induced velocities (perturbations or total,
    ! depending on method) based on given five dipole values (gama), their
    ! locations (x,y,z) and local unit vectors (l and p).
    ! Method is based on interpolating 5 points into linear surface
    ! with equation gama=A*x + B*y + C.
    ! Induced velocities in l and p directions will be coefficients A and B.
    
    subroutine linear_5( ql,qp                                              ,&
                       & x,y,z                                              ,&
                       & l,p                                                ,&
                       & n                                                  ,&
                       & v_err                                              ,&
                       & coplanarity_angle                                  ,&
                       & gama                                               )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                       intent(in)  :: coplanarity_angle
    real(kind=kind_float),    dimension(5),      intent(in)  :: x,y,z,gama
    real(kind=kind_float),    dimension(3),      intent(in)  :: l,p
    real(kind=kind_float),    dimension(3,5),    intent(in)  :: n
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                       intent(out) :: ql,qp
    ! INTENTS INOUT ================================================================
    integer,                    intent(out)                  :: v_err
    ! PRIVATE ======================================================================
    integer                                                  :: i,info,lwork
    real(kind=kind_float)                                    :: correcting_factor
    real(kind=kind_float),    dimension(3)                   :: radius_vector
    real(kind=kind_float),    dimension(6)                   :: work
    real(kind=kind_float),    dimension(5,3)                 :: T
    
    parameter (lwork=6)
    ! ==============================================================================

    ! finding point coordinates in local coordinate system and corresponding
    ! "T" coefficients.
    T(1,:)=(/0._kind_float, 0._kind_float, 1._kind_float/)
    do i=2,5    ! for each of five gama points
        ! find radius vector of "i" point relative to first "1" which is origin of
        ! local coordinate system
        radius_vector=(/x(i)-x(1), y(i)-y(1), z(i)-z(1)/)
        
        ! calculate correction factor that will take into account geometry curvature
        call curvature_correction( correcting_factor                        ,&
                                 & coplanarity_angle                        ,&
                                 & radius_vector                            ,&
                                 & n(:,1)                                   ,&
                                 & n(:,i)                                   )
        
        ! transforming this vector to local coordinate system
        T(i,1)=sum(radius_vector*l)         ! "l" coordinate in local coordinate system (coeff A)
        T(i,2)=sum(radius_vector*p)         ! "p" coordinate in local coordinate system (coeff B)
        T(i,3)=1._kind_float              ! (coeff C)
    enddo

    ! solving system of equations and returning coefficients in gama vector
    if (kind_float .eq. kind_float_l) then
        call sgels( "N"                                                         ,&
                  & 5                                                           ,&
                  & 3                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 5                                                           ,&
                  & gama                                                        ,&
                  & 5                                                           ,&
                  & work                                                        ,&
                  & lwork                                                       ,&
                  & info                                                        )
    else
        call dgels( "N"                                                         ,&
                  & 5                                                           ,&
                  & 3                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 5                                                           ,&
                  & gama                                                        ,&
                  & 5                                                           ,&
                  & work                                                        ,&
                  & lwork                                                       ,&
                  & info                                                        )
    end if    

    if (info .eq. 0) then
        ! localy induced velocities ql,qp
        ql=gama(1)    ! A coefficient
        qp=gama(2)    ! B coefficient
        v_err=0
    else
        ! problem with interpolation occured
        ql=0._kind_float
        qp=0._kind_float
        v_err=1
    endif

    return
    end subroutine linear_5

! SUBROUTINE LINEAR_4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! This function calculates induced velocities (perturbations or total,
    ! depending on method) based on given four dipole values (gama), their
    ! locations (x,y,z) and local unit vectors (l and p).
    ! Method is based on interpolating 4 points into linear surface
    ! with equation gama=A*x + B*y + C.
    ! Induced velocities in l and p directions will be coefficients A and B.
    
    subroutine linear_4( ql,qp                                              ,&
                       & x,y,z                                              ,&
                       & l,p                                                ,&
                       & n                                                  ,&
                       & v_err                                              ,&
                       & coplanarity_angle                                  ,&
                       & gama                                               )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                       intent(in)  :: coplanarity_angle
    real(kind=kind_float),    dimension(4),      intent(in)  :: x,y,z,gama
    real(kind=kind_float),    dimension(3),      intent(in)  :: l,p
    real(kind=kind_float),    dimension(3,4),    intent(in)  :: n
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                       intent(out) :: ql,qp
    ! INTENTS INOUT ================================================================
    integer,                    intent(out)                  :: v_err
    ! PRIVATE ======================================================================
    integer                                                  :: i,info,lwork
    real(kind=kind_float)                                    :: correcting_factor
    real(kind=kind_float),    dimension(3)                   :: radius_vector
    real(kind=kind_float),    dimension(6)                   :: work
    real(kind=kind_float),    dimension(4,3)                 :: T
    
    parameter (lwork=6)
    ! ==============================================================================

    ! finding point coordinates in local coordinate system and corresponding
    ! "T" coefficients.
    T(1,:)=(/0._kind_float, 0._kind_float, 1._kind_float/)
    do i=2,4    ! for each of four gama points
        ! find radius vector of "i" point relative to first "1" which is origin of
        ! local coordinate system
        radius_vector=(/x(i)-x(1), y(i)-y(1), z(i)-z(1)/)
        
        ! calculate correction factor that will take into account geometry curvature
        call curvature_correction( correcting_factor                        ,&
                                 & coplanarity_angle                        ,&
                                 & radius_vector                            ,&
                                 & n(:,1)                                   ,&
                                 & n(:,i)                                   )
        
        ! transforming this vector to local coordinate system
        T(i,1)=sum(radius_vector*l)         ! "l" coordinate in local coordinate system (coeff A)
        T(i,2)=sum(radius_vector*p)         ! "p" coordinate in local coordinate system (coeff B)
        T(i,3)=1._kind_float              ! (coeff C)
    enddo

    ! solving system of equations and returning coefficients in gama vector
    if (kind_float .eq. kind_float_l) then
        call sgels( "N"                                                         ,&
                  & 4                                                           ,&
                  & 3                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 4                                                           ,&
                  & gama                                                        ,&
                  & 4                                                           ,&
                  & work                                                        ,&
                  & lwork                                                       ,&
                  & info                                                        )
    else
        call dgels( "N"                                                         ,&
                  & 4                                                           ,&
                  & 3                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 4                                                           ,&
                  & gama                                                        ,&
                  & 4                                                           ,&
                  & work                                                        ,&
                  & lwork                                                       ,&
                  & info                                                        )
    end if

    if (info .eq. 0) then
        ! localy induced velocities ql,qp
        ql=gama(1)    ! A coefficient
        qp=gama(2)    ! B coefficient
        v_err=0
    else
        ! problem with interpolation occured
        ql=0._kind_float
        qp=0._kind_float
        v_err=1
    endif

    return
    end subroutine linear_4

! SUBROUTINE LINEAR_3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! This function calculates induced velocities (perturbations or total,
    ! depending on method) based on given three dipole values (gama), their
    ! locations (x,y,z) and local unit vectors (l and p).
    ! Method is based on interpolating 3 points into linear surface
    ! with equation gama=A*x + B*y + C.
    ! Induced velocities in l and p directions will be coefficients A and B.
    
    subroutine linear_3( ql,qp                                              ,&
                       & x,y,z                                              ,&
                       & l,p                                                ,&
                       & n                                                  ,&
                       & v_err                                              ,&
                       & coplanarity_angle                                  ,&
                       & gama                                               )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                       intent(in)  :: coplanarity_angle
    real(kind=kind_float),    dimension(3),      intent(in)  :: x,y,z,gama
    real(kind=kind_float),    dimension(3),      intent(in)  :: l,p
    real(kind=kind_float),    dimension(3,3),    intent(in)  :: n
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                       intent(out) :: ql,qp
    ! INTENTS INOUT ================================================================
    integer,                    intent(out)                  :: v_err
    ! PRIVATE ======================================================================
    integer                                                  :: i,info
    integer, dimension(3,3)                                  :: ipiv
    real(kind=kind_float)                                    :: correcting_factor
    real(kind=kind_float),    dimension(3)                   :: radius_vector
    real(kind=kind_float),    dimension(3,3)                 :: T
    ! ==============================================================================

    ! finding point coordinates in local coordinate system and corresponding
    ! "T" coefficients.
    T(1,:)=(/0._kind_float, 0._kind_float, 1._kind_float/)
    do i=2,3    ! for each of three gama points
        ! find radius vector of "i" point relative to first "1" which is origin of
        ! local coordinate system
        radius_vector=(/x(i)-x(1), y(i)-y(1), z(i)-z(1)/)
        
        ! calculate correction factor that will take into account geometry curvature
        call curvature_correction( correcting_factor                        ,&
                                 & coplanarity_angle                        ,&
                                 & radius_vector                            ,&
                                 & n(:,1)                                   ,&
                                 & n(:,i)                                   )
        
        ! transforming this vector to local coordinate system
        T(i,1)=sum(radius_vector*l)         ! "l" coordinate in local coordinate system (coeff A)
        T(i,2)=sum(radius_vector*p)         ! "p" coordinate in local coordinate system (coeff B)
        T(i,3)=1._kind_float              ! (coeff C)
    enddo

    ! solving system of equations and returning coefficients in gama vector
    if (kind_float .eq. kind_float_l) then
        call sgesv( 3                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 3                                                           ,&
                  & ipiv                                                        ,&
                  & gama                                                        ,&
                  & 3                                                           ,&
                  & info                                                        )
    else
        call dgesv( 3                                                           ,&
                  & 1                                                           ,&
                  & T                                                           ,&
                  & 3                                                           ,&
                  & ipiv                                                        ,&
                  & gama                                                        ,&
                  & 3                                                           ,&
                  & info                                                        )
    end if

    if (info .eq. 0) then
        ! localy induced velocities ql,qp
        ql=gama(1)    ! A coefficient
        qp=gama(2)    ! B coefficient
        v_err=0
    else
        ! problem with interpolation occured
        ql=0._kind_float
        qp=0._kind_float
        v_err=1
    endif
    
    return
    end subroutine linear_3
    
! SUBROUTINE curvature_correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! This function calculates correction factor that will take into account geometry 
    ! curvature. Procedure is to assemble a triangle using line that connects center
    ! points (radius_vector) and two angles based on panel normals (normal_1 and normal_2)
    
    subroutine curvature_correction( correcting_factor                      ,&
                                   & coplanarity_angle                      ,&
                                   & radius_vector                          ,&
                                   & normal_1                               ,&
                                   & normal_2                               )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                  intent(in)  :: coplanarity_angle
    real(kind=kind_float),    dimension(3), intent(in)  :: radius_vector                     ,&
                                                         & normal_1                          ,&
                                                         & normal_2
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                  intent(out) :: correcting_factor
    ! PRIVATE ======================================================================
    real(kind=kind_float)                               :: side_c                            ,&
                                                         & average_normal_mod                ,&
                                                         & scalar                            ,&
                                                         & ang1,ang2,ang3
    real(kind=kind_float),    dimension(3)              :: radius_unit_vector                ,&
                                                         & average_normal                    ,&
                                                         & plane_unit_vector                 ,&
                                                         & normal_1_proj                     ,&
                                                         & normal_2_proj
    ! ==============================================================================

    ! calculating radius vector length (length of triangle base line)
    call func_vect_mod(side_c, radius_vector)
    
    ! calculating unit vector of radius_vector
    radius_unit_vector=radius_vector/side_c
    
    ! calculating average normal that will be used for defining projecting plane
    average_normal=(normal_1+normal_2)/2
    ! calculating modulus of that vector
    call func_vect_mod(average_normal_mod, average_normal)
    ! normalizing to unit vector
    average_normal=average_normal/average_normal_mod
    
    ! calculate projection plane (normal unit vector) that will contain euclidian triangle
    call func_vect_prod(plane_unit_vector, radius_unit_vector, average_normal)
    
    ! projecting normals to projection plane
    call func_vect_inprod(scalar,normal_2,plane_unit_vector)
    normal_2_proj=normal_2-scalar*plane_unit_vector
    call func_vect_inprod(scalar,normal_1,plane_unit_vector)
    normal_1_proj=normal_1-scalar*plane_unit_vector
    
    ! calculating first angle of triangle
    call func_vect_ang(radius_unit_vector, normal_2_proj, ang1)
    ang1=abs(ang1-pi/2)
    ! calculating second angle of triangle
    call func_vect_ang(radius_unit_vector, normal_1_proj, ang2)
    ang2=abs(ang2-pi/2)
    ! calculating third angle of triangle
    ang3=pi-ang1-ang2
    
    if (ang3 .lt. coplanarity_angle .and. ang3 .gt. pi-coplanarity_angle) then
        ! calculate triangle sides using sine law
        correcting_factor=abs((sin(ang1)+sin(ang2))/sin(ang3))
    else
        correcting_factor=1.0_kind_float
    endif
    
    return
    end subroutine curvature_correction
    
end module module_velocity
