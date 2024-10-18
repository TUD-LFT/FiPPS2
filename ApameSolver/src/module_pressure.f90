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

! file module_pressure.f90

! This module calculates pressure distributions and sums forces

module module_pressure

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                                            :: press_err
real(kind=kind_float), allocatable, dimension(:,:) :: cp                                     ,&
                                                    & p_dyna                                 ,&
                                                    & p_mano                                 ,&
                                                    & p_stat                                 ,&
                                                    & cl_strip                               ,&
                                                    & cllength_strip                         ,&
                                                    & strip_center
real(kind=kind_float), allocatable, dimension(:)   :: Ftot                                   ,&
                                                    & Fx,Fy,Fz                               ,&
                                                    & Fl,Fm,Fn                               ,&
                                                    & Fdrag,Fside,Flift                      ,&
                                                    & coef_x,coef_y,coef_z                   ,&
                                                    & coef_l,coef_m,coef_n                   ,&
                                                    & coef_drag,coef_side,coef_lift          ,&
                                                    & cl_mom_corrected
integer              , allocatable, dimension(:)   :: elem_strip_id
! ==============================================================================

contains

! SUBROUTINE PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine pressure( interactive                                        ,&
                       & calc_cl_strip                                      ,&
                       & symmetry                                           ,&
                       & node_num                                           ,&
                       & panel_num                                          ,&
                       & panel_num_no_wake                                  ,&
                       & panel_type                                         ,&
                       & case_num                                           ,&
                       & alfa,beta                                          ,&
                       & mac                                                ,&
                       & wing_span                                          ,&
                       & wing_surf                                          ,&
                       & origin                                             ,&
                       & speed                                              ,&
                       & mach                                               ,&
                       & ro                                                 ,&
                       & p_ref                                              ,&
                       & S                                                  ,&
                       & v                                                  ,&
                       & n1,n2,n3                                           ,&
                       & cx,cy,cz                                           ,&
                       & x,y,z                                              ,&
                       & node1,node2,node3,node4                            ,&
                       & elem1,elem2,elem3,elem4                            ,&
                       & error                                              ,&
                       & r_stabilizer                                        )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                                intent(in)                  :: node_num          ,&
                                                                         & panel_num         ,&
                                                                         & panel_num_no_wake ,&
                                                                         & interactive       ,&
                                                                         & calc_cl_strip     ,&
                                                                         & symmetry          ,&
                                                                         & case_num
    integer, dimension(panel_num),          intent(in)                  :: panel_type
    real(kind=kind_float),                                   intent(in) :: mac               ,&
                                                                         & wing_span         ,&
                                                                         & wing_surf         ,&
                                                                         & speed             ,&
                                                                         & mach              ,&
                                                                         & ro                ,&
                                                                         & p_ref             ,&
                                                                         & error             ,&
                                                                         & r_stabilizer
    real(kind=kind_float),    dimension(panel_num),          intent(in) :: S                 ,&
                                                                         & n1,n2,n3          ,&
                                                                         & cx,cy,cz
    real(kind=kind_float),    dimension(node_num),           intent(in) :: x,y,z
    integer, dimension(panel_num),          intent(in)                  :: node1,node2       ,&
                                                                         & node3,node4
    integer, dimension(panel_num),          intent(inout)               :: elem1,elem2       ,&
                                                                         & elem3,elem4
    real(kind=kind_float),    dimension(case_num),           intent(in) :: alfa,beta
    real(kind=kind_float),    dimension(3),                  intent(in) :: origin
    real(kind=kind_float),    dimension(panel_num_no_wake,&
                      &case_num),           intent(in)                  :: v
    ! PRIVATE ======================================================================
    integer                                                             :: i,j,k             ,&
                                                                         & alloc_stat
    real(kind=kind_float)                                               :: q                 ,&
                                                                         & dx,dy,dz          ,&
                                                                         & v_square          ,&
                                                                         & speed_square      ,&
                                                                         & beta_mach         ,&
                                                                         & strip_surface     ,&
                                                                         & strip_length
    integer                                                             :: currEID           ,&
                                                                         & prevEID           ,&
                                                                         & strip_end         ,&
                                                                         & next_panel_case
    integer              ,    dimension(4)                              :: neighbor_panels
    real(kind=kind_float),    dimension(4,3)                            :: strip_nodes
    real(kind=kind_float),    dimension(3)                              :: vec_normal        ,&
                                                                         & Fstrip
    real(kind=kind_float),    dimension(2)                              :: distance
    real(kind=kind_float)                                               :: vec1_mod          ,&
                                                                         & vec2_mod
    ! ==============================================================================
    
    speed_square=speed**2
    q=0.5_kind_float*ro*speed_square
    if (mach .gt. 0._kind_float) then
        beta_mach=sqrt(abs(1._kind_float-mach**2))
    endif
    
    allocate (     cp(panel_num_no_wake,case_num)             ,&
             & p_dyna(panel_num_no_wake,case_num)             ,&
             & p_mano(panel_num_no_wake,case_num)             ,&
             & p_stat(panel_num_no_wake,case_num)             ,&
             & Ftot(case_num)                                 ,&
             & Fx(case_num)                                   ,&
             & Fy(case_num)                                   ,&
             & Fz(case_num)                                   ,&
             & Fl(case_num)                                   ,&
             & Fm(case_num)                                   ,&
             & Fn(case_num)                                   ,&
             & Fdrag(case_num)                                ,&
             & Fside(case_num)                                ,&
             & Flift(case_num)                                ,&
             & coef_x(case_num)                               ,&
             & coef_y(case_num)                               ,&
             & coef_z(case_num)                               ,&
             & coef_l(case_num)                               ,&
             & coef_m(case_num)                               ,&
             & coef_n(case_num)                               ,&
             & coef_drag(case_num)                            ,&
             & coef_side(case_num)                            ,&
             & coef_lift(case_num)                            ,&
             & cl_mom_corrected(case_num)                     ,&
             & stat=alloc_stat                                 )
    if (alloc_stat .ne. 0) then
        call func_message( interactive                                          ,&
                         & "    ERROR: Not enough memory for allocating pressure and force arrays")
        press_err=1
        goto 999
    else
        press_err=0
    endif

    if (calc_cl_strip .eq. 1) then
        allocate ( cl_strip(panel_num-panel_num_no_wake,case_num)           ,&
                 & cllength_strip(panel_num-panel_num_no_wake,case_num)     ,&
                 & strip_center(panel_num-panel_num_no_wake,3)              ,&
                 & elem_strip_id(panel_num_no_wake)                         ,&
                 & stat=alloc_stat                                           )
    else
        allocate ( cl_strip(1,1)                                            ,&
                 & cllength_strip(1,1)                                      ,&
                 & strip_center(1,1)                                        ,&
                 & elem_strip_id(1)                                         ,&
                 & stat=alloc_stat                                           )
    end if
    if (alloc_stat .ne. 0) then
       call func_message( interactive                                          ,&
                        & "    ERROR: Not enough memory for allocating pressure and force arrays")
        press_err=1
        goto 999
    end if
    
    do j=1,case_num
        Fx(j)=0._kind_float
        Fy(j)=0._kind_float
        Fz(j)=0._kind_float
        Fl(j)=0._kind_float
        Fm(j)=0._kind_float
        Fn(j)=0._kind_float
        Ftot(j)=0._kind_float
        do i=1,panel_num_no_wake
            ! avoid dummy panels
            if (panel_type(i) .ne. 20 .and. panel_type(i) .ne. 21) then
                v_square=v(i,j)**2
                
                ! calculating pressure coefficient
                if (mach .gt. 0._kind_float) then
                    cp(i,j)=(1._kind_float-v_square/speed_square)/beta_mach
                else
                    cp(i,j)=1._kind_float-v_square/speed_square
                endif
                
                ! calculating pressures
                p_dyna(i,j)=0.5_kind_float*ro*v_square
                p_mano(i,j)=cp(i,j)*q
                p_stat(i,j)=p_mano(i,j)+p_ref
                
                ! calculating force components
                dx=-p_mano(i,j)*S(i)*n1(i)
                dy=-p_mano(i,j)*S(i)*n2(i)
                dz=-p_mano(i,j)*S(i)*n3(i)
                
                ! summing force components
                Ftot(j)=Ftot(j) + abs(p_mano(i,j))*S(i)
                
                Fx(j)=Fx(j) + dx
                Fy(j)=Fy(j) + dy
                Fz(j)=Fz(j) + dz

                Fl(j)=Fl(j)                        - dy*(cz(i)-origin(3)) + dz*(cy(i)-origin(2))
                Fm(j)=Fm(j) + dx*(cz(i)-origin(3))                        - dz*(cx(i)-origin(1))
                Fn(j)=Fn(j) - dx*(cy(i)-origin(2)) + dy*(cx(i)-origin(1))
                
                ! if symmetry about x-z-plane exists
                if (symmetry .eq. 1) then
                    Fx(j)=Fx(j) + dx
                    Fy(j)=Fy(j) - dy
                    Fz(j)=Fz(j) + dz

                    Fl(j)=Fl(j)                         + dy*(cz(i)-origin(3)) + dz*(-cy(i)-origin(2))
                    Fm(j)=Fm(j) + dx*( cz(i)-origin(3))                        - dz*( cx(i)-origin(1))
                    Fn(j)=Fn(j) - dx*(-cy(i)-origin(2)) - dy*(cx(i)-origin(1))
                end if
            endif
        enddo
        
        ! calculating coefficients in body reference frame
        coef_x(j)=Fx(j)/q/wing_surf
        coef_y(j)=Fy(j)/q/wing_surf
        coef_z(j)=Fz(j)/q/wing_surf
        coef_l(j)=Fl(j)/q/wing_surf/wing_span
        coef_m(j)=Fm(j)/q/wing_surf/mac
        coef_n(j)=Fn(j)/q/wing_surf/wing_span
        
        ! calculating forces and coefficients in aerodynamic reference frame
        Fdrag(j)=  Fx(j)*cos(alfa(j))*cos(beta(j)) + Fy(j)*sin(beta(j)) + Fz(j)*sin(alfa(j))*cos(beta(j))
        Fside(j)=  Fx(j)*cos(alfa(j))*sin(beta(j)) + Fy(j)*cos(beta(j)) + Fz(j)*sin(alfa(j))*sin(beta(j))
        Flift(j)=- Fx(j)*sin(alfa(j))                                   + Fz(j)*cos(alfa(j))
        coef_drag(j)=Fdrag(j)/q/wing_surf
        coef_side(j)=Fside(j)/q/wing_surf
        coef_lift(j)=Flift(j)/q/wing_surf
        
        ! calculating lift coefficient under consideration of the horizontal stabilizer's lift
        ! necessary to compensate the resulting aerodynamic moment
        if (abs(r_stabilizer) .gt. error) then
            cl_mom_corrected(j) = (Flift(j) + Fm(j)/r_stabilizer)/q/wing_surf
        else
            cl_mom_corrected(j) = coef_lift(j)
        end if

        ! calculate stripwise lift coefficient over spanwidth, if requested
        if (calc_cl_strip .eq. 1) then
            cl_strip(:,j) = 0.d0
            elem_strip_id(:) = 0
            ! loop over wake panels -> define the spanwise strips
            do i=panel_num_no_wake+1,panel_num
                strip_end   = elem2(i)
                prevEID     = i
                currEID     = elem1(i)
                strip_nodes(1,:) = (/x(node1(i)), y(node1(i)), z(node1(i))/)
                strip_nodes(2,:) = (/x(node2(i)), y(node2(i)), z(node2(i))/)
                distance(:) = 0.d0
                Fstrip(:) = 0.d0
                ! go stripwise over non-wake panels
                do while (.true.)
                    dx=-p_mano(currEID,j)*S(currEID)*n1(currEID)
                    dy=-p_mano(currEID,j)*S(currEID)*n2(currEID)
                    dz=-p_mano(currEID,j)*S(currEID)*n3(currEID)
                    Fstrip(1)=Fstrip(1) + dx
                    Fstrip(2)=Fstrip(2) - dy
                    Fstrip(3)=Fstrip(3) + dz
                    ! save number of strip for element
                    elem_strip_id(currEID) = i-panel_num_no_wake
                    ! check, if end of strip is reached
                    if (currEID == strip_end) then
                        exit
                    end if
                    ! find out, on which element edge the previous panel is connected
                    if (((node1(currEID) == node1(prevEID)) .or. (node1(currEID) == node2(prevEID)) .or. (node1(currEID) == node3(prevEID)) .or. (node1(currEID) == node4(prevEID))) .and. ((node2(currEID) == node1(prevEID)) .or. (node2(currEID) == node2(prevEID)) .or. (node2(currEID) == node3(prevEID)) .or. (node2(currEID) == node4(prevEID)))) then
                        next_panel_case = 1
                    else if (((node2(currEID) == node1(prevEID)) .or. (node2(currEID) == node2(prevEID)) .or. (node2(currEID) == node3(prevEID)) .or. (node2(currEID) == node4(prevEID))) .and. ((node3(currEID) == node1(prevEID)) .or. (node3(currEID) == node2(prevEID)) .or. (node3(currEID) == node3(prevEID)) .or. (node3(currEID) == node4(prevEID)))) then
                        next_panel_case = 2
                    else if (((node3(currEID) == node1(prevEID)) .or. (node3(currEID) == node2(prevEID)) .or. (node3(currEID) == node3(prevEID)) .or. (node3(currEID) == node4(prevEID))) .and. ((node4(currEID) == node1(prevEID)) .or. (node4(currEID) == node2(prevEID)) .or. (node4(currEID) == node3(prevEID)) .or. (node4(currEID) == node4(prevEID)))) then
                        next_panel_case = 3
                    else if (((node4(currEID) == node1(prevEID)) .or. (node4(currEID) == node2(prevEID)) .or. (node4(currEID) == node3(prevEID)) .or. (node4(currEID) == node4(prevEID))) .and. ((node1(currEID) == node1(prevEID)) .or. (node1(currEID) == node2(prevEID)) .or. (node1(currEID) == node3(prevEID)) .or. (node1(currEID) == node4(prevEID)))) then
                        next_panel_case = 4
                    end if
                    ! look for next panel out of neighboring panels
                    neighbor_panels = (/elem1(currEID),elem2(currEID),elem3(currEID),elem4(currEID)/)
                    do k=1,4
                        if (neighbor_panels(k) .eq. 0) then
                            continue
                        end if                
                        if (next_panel_case == 1) then
                            if (((node3(currEID) == node1(neighbor_panels(k))) .or. (node3(currEID) == node2(neighbor_panels(k))) .or. (node3(currEID) == node3(neighbor_panels(k))) .or. (node3(currEID) == node4(neighbor_panels(k)))) .and. ((node4(currEID) == node1(neighbor_panels(k))) .or. (node4(currEID) == node2(neighbor_panels(k))) .or. (node4(currEID) == node3(neighbor_panels(k))) .or. (node4(currEID) == node4(neighbor_panels(k))))) then
                                prevEID = currEID
                                currEID = neighbor_panels(k)
                                ! check if nodes between current and next panel define the strip corners
                                call func_vect_mod(vec1_mod, (/x(node4(currEID))-strip_nodes(1,1), y(node4(currEID))-strip_nodes(1,2), z(node4(currEID))-strip_nodes(1,3)/))
                                if (vec1_mod .gt. distance(1)) then
                                    distance(1) = vec1_mod
                                    strip_nodes(3,:) = (/x(node4(currEID)), y(node4(currEID)), z(node4(currEID))/)
                                end if
                                call func_vect_mod(vec1_mod, (/x(node3(currEID))-strip_nodes(2,1), y(node3(currEID))-strip_nodes(2,2), z(node3(currEID))-strip_nodes(2,3)/))
                                if (vec1_mod .gt. distance(2)) then
                                    distance(2) = vec1_mod
                                    strip_nodes(4,:) = (/x(node3(currEID)), y(node3(currEID)), z(node3(currEID))/)
                                end if
                                exit
                            end if
                        else if (next_panel_case == 2) then
                            if (((node4(currEID) == node1(prevEID)) .or. (node4(currEID) == node2(prevEID)) .or. (node4(currEID) == node3(prevEID)) .or. (node4(currEID) == node4(prevEID))) .and. ((node1(currEID) == node1(prevEID)) .or. (node1(currEID) == node2(prevEID)) .or. (node1(currEID) == node3(prevEID)) .or. (node1(currEID) == node4(prevEID)))) then
                                prevEID = currEID
                                currEID = neighbor_panels(k)
                                ! check if nodes between current and next panel define the strip corners
                                call func_vect_mod(vec1_mod, (/x(node1(currEID))-strip_nodes(1,1), y(node1(currEID))-strip_nodes(1,2), z(node1(currEID))-strip_nodes(1,3)/))
                                if (vec1_mod .gt. distance(1)) then
                                    distance(1) = vec1_mod
                                    strip_nodes(3,:) = (/x(node1(currEID)), y(node1(currEID)), z(node1(currEID))/)
                                end if
                                call func_vect_mod(vec1_mod, (/x(node4(currEID))-strip_nodes(2,1), y(node4(currEID))-strip_nodes(2,2), z(node4(currEID))-strip_nodes(2,3)/))
                                if (vec1_mod .gt. distance(2)) then
                                    distance(2) = vec1_mod
                                    strip_nodes(4,:) = (/x(node4(currEID)), y(node4(currEID)), z(node4(currEID))/)
                                end if
                                exit
                            end if
                        else if (next_panel_case == 3) then
                            if (((node1(currEID) == node1(prevEID)) .or. (node1(currEID) == node2(prevEID)) .or. (node1(currEID) == node3(prevEID)) .or. (node1(currEID) == node4(prevEID))) .and. ((node2(currEID) == node1(prevEID)) .or. (node2(currEID) == node2(prevEID)) .or. (node2(currEID) == node3(prevEID)) .or. (node2(currEID) == node4(prevEID)))) then
                                prevEID = currEID
                                currEID = neighbor_panels(k)
                                ! check if nodes between current and next panel define the strip corners
                                call func_vect_mod(vec1_mod, (/x(node2(currEID))-strip_nodes(1,1), y(node2(currEID))-strip_nodes(1,2), z(node2(currEID))-strip_nodes(1,3)/))
                                if (vec1_mod .gt. distance(1)) then
                                    distance(1) = vec1_mod
                                    strip_nodes(3,:) = (/x(node2(currEID)), y(node2(currEID)), z(node2(currEID))/)
                                end if
                                call func_vect_mod(vec1_mod, (/x(node1(currEID))-strip_nodes(2,1), y(node1(currEID))-strip_nodes(2,2), z(node1(currEID))-strip_nodes(2,3)/))
                                if (vec1_mod .gt. distance(2)) then
                                    distance(2) = vec1_mod
                                    strip_nodes(4,:) = (/x(node1(currEID)), y(node1(currEID)), z(node1(currEID))/)
                                end if
                                exit
                            end if
                        else if (next_panel_case == 4) then
                            if (((node2(currEID) == node1(prevEID)) .or. (node2(currEID) == node2(prevEID)) .or. (node2(currEID) == node3(prevEID)) .or. (node2(currEID) == node4(prevEID))) .and. ((node3(currEID) == node1(prevEID)) .or. (node3(currEID) == node2(prevEID)) .or. (node3(currEID) == node3(prevEID)) .or. (node3(currEID) == node4(prevEID)))) then
                                prevEID = currEID
                                currEID = neighbor_panels(k)
                                ! check if nodes between current and next panel define the strip corners
                                call func_vect_mod(vec1_mod, (/x(node3(currEID))-strip_nodes(1,1), y(node3(currEID))-strip_nodes(1,2), z(node3(currEID))-strip_nodes(1,3)/))
                                if (vec1_mod .gt. distance(1)) then
                                    distance(1) = vec1_mod
                                    strip_nodes(3,:) = (/x(node3(currEID)), y(node3(currEID)), z(node3(currEID))/)
                                end if
                                call func_vect_mod(vec1_mod, (/x(node2(currEID))-strip_nodes(2,1), y(node2(currEID))-strip_nodes(2,2), z(node2(currEID))-strip_nodes(2,3)/))
                                if (vec1_mod .gt. distance(2)) then
                                    distance(2) = vec1_mod
                                    strip_nodes(4,:) = (/x(node2(currEID)), y(node2(currEID)), z(node2(currEID))/)
                                end if
                                exit
                            end if
                        end if
                    end do
                end do         
                ! calculating vector normal to diagonals of strip (vector product)
                call func_vect_prod(vec_normal,strip_nodes(3,:)-strip_nodes(1,:),strip_nodes(4,:)-strip_nodes(2,:))
                ! calculating modulus of that vector
                call func_vect_mod(strip_surface,vec_normal)
                ! calculated projected area of strip
                strip_surface = strip_surface/2.d0
                ! calculate chord length of stip
                call func_vect_mod(vec1_mod,strip_nodes(4,:)-strip_nodes(1,:))
                call func_vect_mod(vec2_mod,strip_nodes(3,:)-strip_nodes(2,:))
                strip_length = (vec1_mod + vec2_mod)/2.d0
                ! calculate lift coefficient of strip
                cl_strip(i-panel_num_no_wake,j) = (-Fstrip(1)*sin(alfa(j)) + Fstrip(3)*cos(alfa(j)))/q/strip_surface
                cllength_strip(i-panel_num_no_wake,j) = cl_strip(i-panel_num_no_wake,j)*strip_length
                ! calculate center of strip
                strip_center(i-panel_num_no_wake,1) = sum(strip_nodes(:,1))/4.d0
                strip_center(i-panel_num_no_wake,2) = sum(strip_nodes(:,2))/4.d0
                strip_center(i-panel_num_no_wake,3) = sum(strip_nodes(:,3))/4.d0
            end do
        end if
    enddo
        
999 continue
    
    return
    end subroutine pressure

end module module_pressure
