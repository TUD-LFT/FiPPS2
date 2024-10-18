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

! file module_induced_drag.f90

! This module calculates the induced drag using Trefftz-plane integration.
! The theoretical background of the drag calculation is mainly based on [1, 2]
! [1]: Leyser, J.: Accurate computation of lift and induced drag at lifting
!      surfaces; 14th Applied Aerodynamics Conference, Fluid Dynamics and Co-
!      located Conferences; New Orleans, LA, USA; AIAA-96-2405-CP; 1996; doi:
!      10.2514/6.1996-2405
! [2]: Katz, J. and Plotkin, A.: Low-Speed Aerodynamics; 2nd ed.; New York:
!      Cambridge University Press; 2001; ISBN: 978-0-521-66219-2
! [3]: Smith, S.C. and Kroo, I.M.: Computation of induced drag for elliptical
!      and crescent-shaped wings; Journal of Aircraft 30(4); 1993; doi:10.2514/
!      3.46365

module module_induced_drag

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                                            :: idrag_err
real(kind=kind_float), allocatable, dimension(:)   :: Fidrag                                 ,&
                                                      coef_idrag
! ==============================================================================

contains

! SUBROUTINE INDUCED DRAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine induced_drag( gama                                             ,&
                           & sigma                                            ,&
                           & method                                           ,&
                           & node_num                                         ,&
                           & panel_num                                        ,&
                           & panel_num_no_wake                                ,&
                           & case_num                                         ,&
                           & panel_type                                       ,&
                           & interactive                                      ,&
                           & symmetry                                         ,&
                           & x,y,z                                            ,&
                           & speed                                            ,&
                           & speed_x,speed_y,speed_z                          ,&
                           & ro                                               ,&
                           & wing_surf                                        ,&
                           & node1,node2,node3,node4                          ,&
                           & elem1,elem2,elem3,elem4                          ,& 
                           & cx,cy,cz                                         ,&
                           & l1,l2,l3                                         ,&
                           & p1,p2,p3                                         ,&
                           & n1,n2,n3                                         ,&
                           & x1,x2,x3,x4                                      ,&
                           & y1,y2,y3,y4                                      ,&
                           & d1,d2,d3,d4                                      ,&
                           & error                                            )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                                intent(in)                  :: method            ,&
                                                                         & node_num          ,&
                                                                         & panel_num         ,&
                                                                         & panel_num_no_wake ,&
                                                                         & interactive       ,&
                                                                         & symmetry          ,&
                                                                         & case_num
    real(kind=kind_float),                                   intent(in) :: speed             ,&
                                                                         & wing_surf         ,&
                                                                         & ro                ,&
                                                                         & error
    integer, dimension(panel_num),          intent(in)                  :: panel_type        ,&
                                                                         & node1,node2       ,&
                                                                         & node3,node4       ,&
                                                                         & elem1,elem2       ,&
                                                                         & elem3,elem4
    real(kind=kind_float),    dimension(panel_num),          intent(in) :: cx,cy,cz          ,&
                                                                         & l1,l2,l3          ,&
                                                                         & p1,p2,p3          ,&
                                                                         & n1,n2,n3          ,&
                                                                         & x1,x2,x3,x4       ,&
                                                                         & y1,y2,y3,y4       ,&
                                                                         & d1,d2,d3,d4
    real(kind=kind_float),    dimension(node_num),           intent(in) :: x,y,z
    real(kind=kind_float),    dimension(case_num),           intent(in) :: speed_x ,&
                                                                         & speed_y ,&
                                                                         & speed_z
    real(kind=kind_float),    dimension(panel_num_no_wake,case_num), intent(in) :: gama      ,&
                                                                         & sigma
    ! PRIVATE ======================================================================
    integer                                                             :: h,i,j,k,m         ,&
                                                                         & interp_error      ,&
                                                                         & alloc_stat        ,&
                                                                         & n_panel_wake      ,&
                                                                         & n_intpnt          ,&
                                                                         & pid               ,&
                                                                         & help
    real(kind=kind_float)                                               :: gama_w1, gama_w2  ,&
                                                                         & q                 ,&
                                                                         & sym_fac
    real(kind=kind_float), dimension(3)                                 :: u_ind             ,&
                                                                         & xyz_1             ,&
                                                                         & xyz_2             ,&
                                                                         & xyz_3             ,&
                                                                         & xyz_4             ,&
                                                                         & u_tmp             ,&
                                                                         & point_l           ,&
                                                                         & y_neg
    real(kind=kind_float), dimension(3,3)                               :: TRMat_gl, TRMat_lg
    real(kind=kind_float),    allocatable, dimension(:)                 :: fac               ,&
                                                                         & dl
    real(kind=kind_float),    allocatable, dimension(:,:)               :: dir
    real(kind=kind_float),    allocatable, dimension(:,:,:)             :: point
    real(kind=kind_float), dimension(1)                                 :: w_gauss1, x_gauss1
    real(kind=kind_float), dimension(2)                                 :: w_gauss2, x_gauss2
    real(kind=kind_float), dimension(3)                                 :: w_gauss3, x_gauss3
    real(kind=kind_float), dimension(4)                                 :: w_gauss4, x_gauss4
    parameter (y_neg=(/1._kind_float,-1._kind_float,1._kind_float/))
    parameter (n_intpnt=3)
    parameter (w_gauss1=(/2._kind_float/))
    parameter (x_gauss1=(/0._kind_float/))
    parameter (w_gauss2=(/-1._kind_float/sqrt(3._kind_float) ,&
                         & 1._kind_float/sqrt(3._kind_float)/))
    parameter (x_gauss2=(/ 1._kind_float ,&
                         & 1._kind_float/))
    parameter (w_gauss3=(/(5._kind_float/9._kind_float) ,&
                         &(8._kind_float/9._kind_float) ,&
                         &(5._kind_float/9._kind_float)/))
    parameter (x_gauss3=(/-sqrt(3._kind_float/5._kind_float) ,&
                         & 0._kind_float                     ,&
                         & sqrt(3._kind_float/5._kind_float)/))
    parameter (w_gauss4=(/ (18._kind_float-sqrt(30._kind_float))/36._kind_float   ,&
                         & (18._kind_float+sqrt(30._kind_float))/36._kind_float   ,&
                         & (18._kind_float+sqrt(30._kind_float))/36._kind_float   ,&
                         & (18._kind_float-sqrt(30._kind_float))/36._kind_float/))
    parameter (x_gauss4=(/-sqrt(3._kind_float/7._kind_float+2._kind_float/7._kind_float*sqrt(6._kind_float/5._kind_float))  ,&
                         &-sqrt(3._kind_float/7._kind_float-2._kind_float/7._kind_float*sqrt(6._kind_float/5._kind_float))  ,&
                         & sqrt(3._kind_float/7._kind_float-2._kind_float/7._kind_float*sqrt(6._kind_float/5._kind_float))  ,&
                         & sqrt(3._kind_float/7._kind_float+2._kind_float/7._kind_float*sqrt(6._kind_float/5._kind_float))/))
    ! ==============================================================================
    
    ! presume no interpolation errors occured
    interp_error=0

    ! presume no induced_drag subroutine allocation errors
    idrag_err=0    
    
    ! number of wake-panels
    n_panel_wake = panel_num - panel_num_no_wake
    
    ! allocate vectors of induced velocities and circulations
    allocate(dir(n_panel_wake,3)              ,&
           & dl(n_panel_wake)                 ,&
           & fac(n_intpnt)                    ,&
           & point(n_panel_wake,n_intpnt,3)   ,&
           & Fidrag(case_num)                 ,&
           & coef_idrag(case_num)             ,&
           & stat=alloc_stat)
    if (alloc_stat .ne. 0) then
        call func_message( interactive                                          ,&
                         & "    ERROR: Not enough memory for allocating induced drag arrays")
        idrag_err=1
        goto 999
    else
        idrag_err=0
    endif

    ! calculate factors for Gaussian quadrature in interval [0,1]
    fac(:) = 0.5_kind_float*x_gauss3 + 0.5_kind_float
    
    ! loop over wake panels
    do i=1,n_panel_wake
        ! panel-id
        pid = panel_num_no_wake+i

        ! calculate direction vector along trailing
        ! edge of current wake panel with respect
        ! to global coordinate system
        dir(i,1) = x(node3(pid)) - x(node4(pid))
        dir(i,2) = y(node3(pid)) - y(node4(pid))
        dir(i,3) = z(node3(pid)) - z(node4(pid))

        ! loop over integration points for current wake panel
        do j=1,n_intpnt
            ! calculate coordinates of integration points
            ! in global coordinate system
            point(i,j,1) = x(node4(pid)) + fac(j)*dir(i,1)
            point(i,j,2) = y(node4(pid)) + fac(j)*dir(i,2)
            point(i,j,3) = z(node4(pid)) + fac(j)*dir(i,3)
        end do

        ! Length of line segment
        call func_vect_mod(dl(i), dir(i,:))

        ! Normalize direction vector to unity
        dir(i,:) = dir(i,:) / dl(i)
    end do
    
    ! initialize force vector
    Fidrag(:) = 0._kind_float
    
    ! loop over cases
    do m=1,case_num

        ! loop over wake panels
        do i = 1, n_panel_wake
            ! panel-id
            pid = panel_num_no_wake+i
            
            ! strength of wake panel
            gama_w1 = gama(elem1(pid),m) - gama(elem2(pid),m)

            ! loop over integration points per wake panel
            do j=1,n_intpnt

                ! initialize vector of induced velocity at current
                ! integration point
                u_ind(:) = 0._kind_float

                ! loop over influencing wake panels
                do k = panel_num_no_wake+1, panel_num

                    ! get nodal coordinates of influencing wake panel
                    ! with respect to local coordinate system of
                    ! influencing wake panel
                    xyz_1(1) = x1(k); xyz_1(2) = y1(k); xyz_1(3) = 0._kind_float
                    xyz_2(1) = x2(k); xyz_2(2) = y2(k); xyz_2(3) = 0._kind_float
                    xyz_3(1) = x3(k); xyz_3(2) = y3(k); xyz_3(3) = 0._kind_float
                    xyz_4(1) = x4(k); xyz_4(2) = y4(k); xyz_4(3) = 0._kind_float

                    ! set up transformation matrix for transformation of
                    ! vectors from global to local coordinate system of
                    ! influencing wake panel
                    TRMat_gl(1,:) = (/l1(k),l2(k),l3(k)/)
                    TRMat_gl(2,:) = (/p1(k),p2(k),p3(k)/)
                    TRMat_gl(3,:) = (/n1(k),n2(k),n3(k)/)
                    
                    ! set up transformation matrix for transformation of
                    ! vectors from local coordinate system of influencing
                    ! wake panel to global coordinate system
                    TRMat_lg(:,:) = transpose(TRMat_gl(:,:))
                
                    ! strength of influencing wake panel
                    if (k .gt. panel_num_no_wake) then
                        gama_w2 = gama(elem1(k),m) - gama(elem2(k),m)
                    else
                        gama_w2 = gama(k,m)
                    end if
                
                    ! if no symmetry about x-z-plane must be considered
                    if (symmetry .eq. 0) then
                        help=1
                    ! if symmetry about x-z-plane must be considered
                    else
                        help=2
                    end if
                    
                    do h=1,help

                        ! transformation of current integration point to
                        ! local coordinate system of influencing panel
                        !   1st: subtract center of element
                        !   2nd: rotation via transformation matrix
                        if (h .eq. 1) then
                            point_l(:) = point(i,j,:) - (/cx(k),cy(k),cz(k)/)
                        else
                            point_l(:) = point(i,j,:)*y_neg(:) - (/cx(k),cy(k),cz(k)/)
                        end if
                        point_l = matmul(TRMat_gl,point_l)
                        
                        ! induced velocity of first trailing vortex
                        ! of influencing wake panel with respect to
                        ! local coordinate system
                        call induced_vorticity_vortexline( gama_w2,       &
                                                         & point_l,       &
                                                         & xyz_4,         &
                                                         & xyz_1,         &
                                                         & u_tmp,         &
                                                         & error)

                        ! transform induced velocity to global coordinate
                        ! system and add to induced velocity at current
                        ! integration point
                        if (h .eq. 1) then
                            u_ind(:) = u_ind(:) + matmul(TRMat_lg,u_tmp)
                        else
                            u_ind(:) = u_ind(:) + matmul(TRMat_lg,u_tmp)*y_neg(:)
                        end if
                        
                        ! induced velocity of second trailing vortex
                        ! of influencing wake panel with respect to
                        ! local coordinate system
                        call induced_vorticity_vortexline( gama_w2,       &
                                                         & point_l,       &
                                                         & xyz_2,         &
                                                         & xyz_3,         &
                                                         & u_tmp,         &
                                                         & error)

                        ! transform induced velocity to global coordinate
                        ! system and add to induced velocity at current
                        ! integration point
                        if (h .eq. 1) then
                            u_ind(:) = u_ind(:) + matmul(TRMat_lg,u_tmp)
                        else
                            u_ind(:) = u_ind(:) + matmul(TRMat_lg,u_tmp)*y_neg(:)
                        end if
                    end do
                end do

                ! numerical integration of force via Gaussian quadrature
                if (symmetry .eq. 1) then
                    sym_fac = 2._kind_float
                else
                    sym_fac = 1._kind_float
                end if

                Fidrag(m) = Fidrag(m) + sym_fac*gama_w1*dot_product(u_ind(:),(/n1(pid),n2(pid),n3(pid)/))*dl(i)*w_gauss3(j)/2._kind_float
            end do
        end do

        ! multiply with air density
        Fidrag(m) = Fidrag(m)*ro

        ! dynamic pressure
        q=0.5_kind_float*ro*speed**2

        ! drag coefficient
        coef_idrag(m) = Fidrag(m)/q/wing_surf
    end do

    ! deallocate no longer used arrays
    deallocate(dir     ,&
             & dl      ,&
             & fac     ,&
             & point)

999 continue
    
    return
    end subroutine induced_drag
    
! SUBROUTINE induced_vorticity_vortexline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    subroutine induced_vorticity_vortexline( gamma,        &
                                          & point,         &
                                          & edge_1,        & 
                                          & edge_2,        &
                                          & u_ind,         &
                                          & error)

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                  intent(in)  :: gamma,  &
                                                         & error
    real(kind=kind_float), dimension(3),    intent(in)  :: point,  &
                                                         & edge_1, &
                                                         & edge_2
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), dimension(3),    intent(out) :: u_ind
    ! PRIVATE ======================================================================
    real(kind=kind_float), dimension(3)                 :: r_0, r_1, r_2 , cross_12
    real(kind=kind_float)                               :: abs_r1, abs_r2, abs_cross12
    real(kind=kind_float)                               :: k_fac
    ! ==============================================================================
    
    ! get vector from begin to end of vortexline
    r_0(:) = edge_2(:) - edge_1(:)
    ! get vector from begin of vortexline to point of interest
    r_1(:) = point(:) - edge_1(:)
    ! get vector from end of vortexline to point of interest
    r_2(:) = point(:) - edge_2(:)

    ! get absolute of r_1
    call func_vect_mod(abs_r1, r_1(:))
    ! get absolute of r_2
    call func_vect_mod(abs_r2, r_2(:))

    ! get cross product of r_1 with r_2
    call func_vect_prod(cross_12(:), r_1(:), r_2(:))
    ! get absolute of cross product of r_1 with r_2
    call func_vect_mod(abs_cross12, cross_12(:))

    ! if point of interest lies on vortex line or its elongation
    if ((abs_r1 .lt. error) .or. (abs_r2 .lt. error) .or. (abs_cross12**2 .lt. error)) then
        ! singular solution if point of interest lies on vortex line        
        u_ind(:) = 0._kind_float

    else
        ! induced vorticity according to biot-savart law 
        k_fac = gamma/(4._kind_float*pi*abs_cross12**2)
        k_fac = k_fac*(dot_product(r_0,r_1)/abs_r1 - dot_product(r_0,r_2)/abs_r2)
        
        u_ind = cross_12(:)*k_fac
    end if
    
    end subroutine induced_vorticity_vortexline

! SUBROUTINE induced_vorticity_doublet_panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    subroutine induced_vorticity_doublet_panel( gamma,        &
                                              & point,        &
                                              & node_1,       & 
                                              & node_2,       &
                                              & node_3,       &
                                              & node_4,       &
                                              & u_ind,        &
                                              & error)

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                  intent(in)  :: gamma,  &
                                                         & error
    real(kind=kind_float), dimension(3),    intent(in)  :: point,  &
                                                         & node_1, &
                                                         & node_2, &
                                                         & node_3, &
                                                         & node_4
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), dimension(3),    intent(out) :: u_ind
    ! PRIVATE ======================================================================
    real(kind=kind_float)                               :: fac,          &
                                                         & ax, ay, az,   &
                                                         & bx, by, bz
    real(kind=kind_float), dimension(4,3)               :: nodes
    real(kind=kind_float), dimension(4)                 :: r
    integer                                             :: i,j
    ! ==============================================================================

    ! get absolute of r_1
    call func_vect_mod(r(1), point(:) - node_1(:))
    ! get absolute of r_2
    call func_vect_mod(r(2), point(:) - node_2(:))
    ! get absolute of r_3
    call func_vect_mod(r(3), point(:) - node_3(:))
    ! get absolute of r_4
    call func_vect_mod(r(4), point(:) - node_4(:))

    nodes(1,:) = node_1(:)
    nodes(2,:) = node_2(:)
    nodes(3,:) = node_3(:)
    nodes(4,:) = node_4(:)
    
    u_ind(:) = 0._kind_float
    do i=1,4
        j = i+1
        if (j .gt. 4) j=1

        ax = point(1) - nodes(i,1)
        ay = point(2) - nodes(i,2)
        az = point(3) - nodes(i,3)

        bx = point(1) - nodes(j,1)
        by = point(2) - nodes(j,2)
        bz = point(3) - nodes(j,3)
        
        fac = r(i)*r(j) + ax*bx + ay*by + az*bz
        if (abs(fac) .lt. error) then
            fac = 0._kind_float
        else
            fac = (r(i)+r(j))/(r(i)*r(j)*fac)
        end if
        
        u_ind(1) = u_ind(1) + (ay*bz - az*by)*fac
        u_ind(2) = u_ind(2) + (az*bx - ax*bz)*fac
        u_ind(3) = u_ind(3) + (ax*by - ay*bx)*fac
    end do
    u_ind(:) = u_ind(:)*gamma/(4._kind_float*pi)
    
    end subroutine induced_vorticity_doublet_panel

! SUBROUTINE velocity_potential_doublet_panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    subroutine velocity_potential_doublet_panel( gamma,        &
                                               & point,        &
                                               & node_1,       & 
                                               & node_2,       &
                                               & node_3,       &
                                               & node_4,       &
                                               & phi,          &
                                               & error)

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                  intent(in)  :: gamma,  &
                                                         & error
    real(kind=kind_float), dimension(3),    intent(in)  :: point,  &
                                                         & node_1, &
                                                         & node_2, &
                                                         & node_3, &
                                                         & node_4
    ! INTENTS OUT ==================================================================
    real(kind=kind_float),                  intent(out) :: phi
    ! PRIVATE ======================================================================
    real(kind=kind_float)                               :: m,            &
                                                         & ei, ej,       &
                                                         & hi, hj
    real(kind=kind_float), dimension(4,3)               :: nodes
    real(kind=kind_float), dimension(4)                 :: r
    integer                                             :: i,j
    ! ==============================================================================

    ! get absolute of r_1
    call func_vect_mod(r(1), point(:) - node_1(:))
    ! get absolute of r_2
    call func_vect_mod(r(2), point(:) - node_2(:))
    ! get absolute of r_3
    call func_vect_mod(r(3), point(:) - node_3(:))
    ! get absolute of r_4
    call func_vect_mod(r(4), point(:) - node_4(:))

    nodes(1,:) = node_1(:)
    nodes(2,:) = node_2(:)
    nodes(3,:) = node_3(:)
    nodes(4,:) = node_4(:)
    
    phi = 0._kind_float
    do i=1,4
        j = i+1
        if (j .gt. 4) j=1

        ei = (point(1) - nodes(i,1))**2.d0 + point(3)**2.d0
        ej = (point(1) - nodes(j,1))**2.d0 + point(3)**2.d0

        m  = (nodes(j,2) - nodes(i,2))/(nodes(j,1) - nodes(i,1))
        
        hi = (point(1) - nodes(i,1))*(point(2) - nodes(i,2))
        hj = (point(1) - nodes(j,1))*(point(2) - nodes(j,2))
        
        phi = phi + atan2((m*ei - hi),(point(3)*r(i))) - atan2((m*ej-hj),(point(3)*r(j)))

    end do
    phi = phi*gamma/(4._kind_float*pi)
    
    end subroutine velocity_potential_doublet_panel

! SUBROUTINE induced_vorticity_strength_panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    subroutine induced_vorticity_strength_panel( sigma,        &
                                               & point,        &
                                               & node_1,       & 
                                               & node_2,       &
                                               & node_3,       &
                                               & node_4,       &
                                               & u_ind,        &
                                               & error)

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float),                  intent(in)  :: sigma,  &
                                                         & error
    real(kind=kind_float), dimension(3),    intent(in)  :: point,  &
                                                         & node_1, &
                                                         & node_2, &
                                                         & node_3, &
                                                         & node_4
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), dimension(3),    intent(out) :: u_ind
    ! PRIVATE ======================================================================
    real(kind=kind_float)                               :: fac,          &
                                                         & s, m,         &
                                                         & ei, ej,       &
                                                         & hi, hj
    real(kind=kind_float), dimension(4,3)               :: nodes
    real(kind=kind_float), dimension(4)                 :: r
    integer                                             :: i,j
    ! ==============================================================================

    ! get absolute of r_1
    call func_vect_mod(r(1), point(:) - node_1(:))
    ! get absolute of r_2
    call func_vect_mod(r(2), point(:) - node_2(:))
    ! get absolute of r_3
    call func_vect_mod(r(3), point(:) - node_3(:))
    ! get absolute of r_4
    call func_vect_mod(r(4), point(:) - node_4(:))

    nodes(1,:) = node_1(:)
    nodes(2,:) = node_2(:)
    nodes(3,:) = node_3(:)
    nodes(4,:) = node_4(:)
    
    u_ind(:) = 0._kind_float
    do i=1,4
        j = i+1
        if (j .gt. 4) j=1

        call func_vect_mod(s, nodes(j,:) - nodes(i,:))

        if (((r(i)+r(j)+s) .gt. error) .and. (s .gt. error)) then
            fac = log((r(i)+r(j)-s)/(r(i)+r(j)+s))/s
        else
            fac = 0._kind_float
        end if
        
        u_ind(1) = u_ind(1) + (nodes(j,2) - nodes(i,2))*fac
        u_ind(2) = u_ind(2) + (nodes(i,1) - nodes(j,1))*fac
        
        ei = (point(1) - nodes(i,1))**2.d0 + point(3)**2.d0
        ej = (point(1) - nodes(j,1))**2.d0 + point(3)**2.d0

        m  = (nodes(j,2) - nodes(i,2))/(nodes(j,1) - nodes(i,1))
        
        hi = (point(1) - nodes(i,1))*(point(2) - nodes(i,2))
        hj = (point(1) - nodes(j,1))*(point(2) - nodes(j,2))
        
        u_ind(3) = u_ind(3) + atan2((m*ei - hi),(point(3)*r(i))) - atan2((m*ej-hj),(point(3)*r(j)))

    end do
    u_ind(:) = u_ind(:)*sigma/(4._kind_float*pi)
    
    end subroutine induced_vorticity_strength_panel

end module module_induced_drag
