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

! file module_grid.f90

! This module contains grid subroutine that calculates grid information based on
! nodes and elements. Grid information consists of:
!   - farfield distances (FF)
!   - panel surfaces (S)
!   - components of unit normal vector (n1,n2,n3)
!   - components of unit lengthwise vector (l1,l2,l3)
!   - components of unit perpendicular vector (p1,p2,p3)
!   - panel center coordinates (cx,cy,cz)
!   - panel collocations point coordinates (colx,coly,colz)
!   - panel node coordinates in local coordinate system (x1,x2,x3,x4,y1,y2,y3,y4)
!   - panel side lengths (d1,d2,d3,d4)

module module_grid

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                                            :: grid_err
real(kind=kind_float),   allocatable, dimension(:) :: FF                                     ,&
                                                    & S                                      ,&
                                                    & n1,n2,n3                               ,&
                                                    & l1,l2,l3                               ,&
                                                    & p1,p2,p3                               ,&
                                                    & cx,cy,cz                               ,&
                                                    & colx,coly,colz                         ,&
                                                    & x1,x2,x3,x4                            ,&
                                                    & y1,y2,y3,y4                            ,&
                                                    & d1,d2,d3,d4
! ==============================================================================

contains

! SUBROUTINE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine grid( interactive                                            ,&
                   & collcalc                                               ,&
                   & panel_num                                              ,&
                   & panel_num_no_wake                                      ,&
                   & node_num                                               ,&
                   & panel_type                                             ,&
                   & node1,node2,node3,node4                                ,&
                   & farfield                                               ,&
                   & error                                                  ,&
                   & colldist                                               ,&
                   & x,y,z                                                  ,&
                   & wake_only                                               )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                       intent(in)                  :: panel_num                  ,&
                                                                & panel_num_no_wake          ,&
                                                                & node_num                   ,&
                                                                & interactive                ,&
                                                                & collcalc
    integer, dimension(panel_num), intent(in)                  :: panel_type                 ,&
                                                                & node1,node2,node3,node4
    real(kind=kind_float),                          intent(in) :: farfield                   ,&
                                                                & error                      ,&
                                                                & colldist
    real(kind=kind_float),    dimension(node_num),  intent(in) :: x,y,z
    logical, intent(in)                                        :: wake_only
    ! PRIVATE ======================================================================
    integer                                   :: i,alloc_stat,ind_start
    ! ==============================================================================
    
    if (wake_only .eqv. .false.) then
        allocate( FF(panel_num)                                                 ,&
                & S(panel_num)                                                  ,&
                & n1(panel_num),n2(panel_num),n3(panel_num)                     ,&
                & l1(panel_num),l2(panel_num),l3(panel_num)                     ,&
                & p1(panel_num),p2(panel_num),p3(panel_num)                     ,&
                & cx(panel_num),cy(panel_num),cz(panel_num)                     ,&
                & x1(panel_num),x2(panel_num),x3(panel_num),x4(panel_num)       ,&
                & y1(panel_num),y2(panel_num),y3(panel_num),y4(panel_num)       ,&
                & d1(panel_num),d2(panel_num),d3(panel_num),d4(panel_num)       ,&
                & colx(panel_num),coly(panel_num),colz(panel_num)               ,&
                & stat=alloc_stat                                               )
        if (alloc_stat .ne. 0) then
            call func_message( interactive                                          ,&
                & "    ERROR: Not enough memory for allocating grid data fields")
            grid_err=1
            goto 999
        endif
    end if
    
    if (wake_only .eq. .true.) then
        ind_start = panel_num_no_wake + 1
    else
        ind_start = 1
    end if
    
    do i=ind_start,panel_num
        ! in case of quadrilateral panel
        if ( panel_type(i) .eq. 1   .or. &
           & panel_type(i) .eq. 10  .or. &
           & panel_type(i) .eq. 20  ) then
            call grid_panel_quad( interactive                                       ,&
                                & collcalc                                          ,&
                                & farfield                                          ,&
                                & error                                             ,&
                                & colldist                                          ,&
                                & x(node1(i)),x(node2(i)),x(node3(i)),x(node4(i))   ,&
                                & y(node1(i)),y(node2(i)),y(node3(i)),y(node4(i))   ,&
                                & z(node1(i)),z(node2(i)),z(node3(i)),z(node4(i))   ,&
                                & FF(i)                                             ,&
                                & S(i)                                              ,&
                                & n1(i),n2(i),n3(i)                                 ,&
                                & l1(i),l2(i),l3(i)                                 ,&
                                & p1(i),p2(i),p3(i)                                 ,&
                                & cx(i),cy(i),cz(i)                                 ,&
                                & colx(i),coly(i),colz(i)                           ,&
                                & x1(i),x2(i),x3(i),x4(i)                           ,&
                                & y1(i),y2(i),y3(i),y4(i)                           ,&
                                & d1(i),d2(i),d3(i),d4(i)                           ,&
                                & grid_err                                          )
        ! in case of triangular panel
        elseif ( panel_type(i) .eq. 2   .or.&
               & panel_type(i) .eq. 11  .or.&
               & panel_type(i) .eq. 21  ) then
            call grid_panel_tri( interactive                                ,&
                               & collcalc                                   ,&
                               & farfield                                   ,&
                               & error                                      ,&
                               & colldist                                   ,&
                               & x(node1(i)),x(node2(i)),x(node3(i))        ,&
                               & y(node1(i)),y(node2(i)),y(node3(i))        ,&
                               & z(node1(i)),z(node2(i)),z(node3(i))        ,&
                               & FF(i)                                      ,&
                               & S(i)                                       ,&
                               & n1(i),n2(i),n3(i)                          ,&
                               & l1(i),l2(i),l3(i)                          ,&
                               & p1(i),p2(i),p3(i)                          ,&
                               & cx(i),cy(i),cz(i)                          ,&
                               & colx(i),coly(i),colz(i)                    ,&
                               & x1(i),x2(i),x3(i)                          ,&
                               & y1(i),y2(i),y3(i)                          ,&
                               & d1(i),d2(i),d3(i)                          ,&
                               & grid_err                                   )
        endif
        if (grid_err .eq. 1) then
            goto 999
        endif
    enddo
    
999 continue
    
    return
    end subroutine grid
    
! SUBROUTINE GRID_PANEL_QUAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! subroutine for calculating grid data of quadrilateral panels
    
    subroutine grid_panel_quad( interactive                                 ,&
                              & collcalc                                    ,&
                              & farfield                                    ,&
                              & error                                       ,&
                              & colldist                                    ,&
                              & x1,x2,x3,x4                                 ,&
                              & y1,y2,y3,y4                                 ,&
                              & z1,z2,z3,z4                                 ,&
                              & FF,S                                        ,&
                              & n1,n2,n3                                    ,&
                              & l1,l2,l3                                    ,&
                              & p1,p2,p3                                    ,&
                              & cx,cy,cz                                    ,&
                              & colx,coly,colz                              ,&
                              & xl1,xl2,xl3,xl4                             ,&
                              & yl1,yl2,yl3,yl4                             ,&
                              & d1,d2,d3,d4                                 ,&
                              & grid_err                                    )

    implicit none
    
    ! INTENTS IN ===================================================================
    integer, intent(in)                    :: interactive, collcalc
    real(kind=kind_float),    intent(in)   :: x1,x2,x3,x4                                    ,&
                                            & y1,y2,y3,y4                                    ,&
                                            & z1,z2,z3,z4                                    ,&
                                            & farfield                                       ,&
                                            & error                                          ,&
                                            & colldist
    ! INTENTS OUT ==================================================================
    integer, intent(out)                   :: grid_err
    real(kind=kind_float),    intent(out)  :: FF                                             ,&
                                            & S                                              ,&
                                            & n1,n2,n3                                       ,&
                                            & l1,l2,l3                                       ,&
                                            & p1,p2,p3                                       ,&
                                            & cx,cy,cz                                       ,&
                                            & colx,coly,colz                                 ,&
                                            & xl1,xl2,xl3,xl4                                ,&
                                            & yl1,yl2,yl3,yl4                                ,&
                                            & d1,d2,d3,d4
    ! PRIVATE ======================================================================
    real(kind=kind_float)                  :: A_mod,B_mod,C_mod                              ,&
                                            & lx,ly,lz                                       ,&
                                            & cx1,cx2,cx3,cx4                                ,&
                                            & cy1,cy2,cy3,cy4                                ,&
                                            & cz1,cz2,cz3,cz4                                ,&
                                            & l_mod                                          ,&
                                            & dsum
    real(kind=kind_float),    dimension(3) :: A,B,C,p
    ! ==============================================================================
    
    grid_err=0

    ! determining diagonals as vectors
    A=(/(x3-x1), (y3-y1), (z3-z1)/)
    B=(/(x4-x2), (y4-y2), (z4-z2)/)
    
    ! calculating their modules
    call func_vect_mod(A_mod,A)
    call func_vect_mod(B_mod,B)
    
    ! calculating far field distance as longer diagonal multiplied by farfield factor
    FF=farfield*max(A_mod,B_mod)
    
    ! calculating vector normal to diagonals (vector product)
    call func_vect_prod(C,A,B)
    
    ! calculating modulus of that vector
    call func_vect_mod(C_mod,C)
    
    ! panel surface is exactly half of that modulus
    S=C_mod/2._kind_float
    if (C_mod .lt. error) then
        call func_message( interactive                                      ,&
                         & "ERROR: Panel with zero or negative surface detected")
        grid_err=1
        goto 888
    endif
    
    ! normalizing normal vector to unit vector
    n1=C(1)/C_mod
    n2=C(2)/C_mod
    n3=C(3)/C_mod
    
    ! calculating lengthwise unit vector (/l1,l2,l3/)
    lx=((x4+x3)-(x1+x2))/2
    ly=((y4+y3)-(y1+y2))/2
    lz=((z4+z3)-(z1+z2))/2
    call func_vect_mod(l_mod,(/lx,ly,lz/))
    if (l_mod .lt. error) then
        call func_message( interactive                                      ,&
                         & "ERROR: Zero length lengthwise unit vector detected")
        grid_err=1
        goto 888
    endif
    l1=lx/l_mod
    l2=ly/l_mod
    l3=lz/l_mod
    
    ! calculating unit vector perpendicular to normal and lengthwise unit vectors (/p1,p2,p3/)
    call func_vect_prod(p,(/n1,n2,n3/),(/l1,l2,l3/))
    p1=p(1)
    p2=p(2)
    p3=p(3)
    
    if (collcalc .eq. 0) then
        ! approximative collocation/center points calculation
        
        ! center points
        cx=(x1 + x2 + x3 + x4)/4
        cy=(y1 + y2 + y3 + y4)/4
        cz=(z1 + z2 + z3 + z4)/4
        
        ! panel node coordinates (in local CS)
        xl1=(x1-cx)*l1 + (y1-cy)*l2 + (z1-cz)*l3
        xl2=(x2-cx)*l1 + (y2-cy)*l2 + (z2-cz)*l3
        xl3=(x3-cx)*l1 + (y3-cy)*l2 + (z3-cz)*l3
        xl4=(x4-cx)*l1 + (y4-cy)*l2 + (z4-cz)*l3
        yl1=(x1-cx)*p1 + (y1-cy)*p2 + (z1-cz)*p3
        yl2=(x2-cx)*p1 + (y2-cy)*p2 + (z2-cz)*p3
        yl3=(x3-cx)*p1 + (y3-cy)*p2 + (z3-cz)*p3
        yl4=(x4-cx)*p1 + (y4-cy)*p2 + (z4-cz)*p3
        
        ! panel side lengths
        d1=sqrt((xl2-xl1)**2 + (yl2-yl1)**2)
        d2=sqrt((xl3-xl2)**2 + (yl3-yl2)**2)
        d3=sqrt((xl4-xl3)**2 + (yl4-yl3)**2)
        d4=sqrt((xl1-xl4)**2 + (yl1-yl4)**2)
    elseif (collcalc .eq. 1) then
        ! accurate collocation/center points calculation
        
        ! panel side lengths and their sum
        call func_vect_mod(d1,(/(x2-x1),(y2-y1),(z2-z1)/))
        call func_vect_mod(d2,(/(x3-x2),(y3-y2),(z3-z2)/))
        call func_vect_mod(d3,(/(x4-x3),(y4-y3),(z4-z3)/))
        call func_vect_mod(d4,(/(x1-x4),(y1-y4),(z1-z4)/))
        dsum=d1 + d2 + d3 + d4
        if (dsum .lt. error) then
            call func_message( interactive                                      ,&
                             & "ERROR: Zero sum of panel sides detected")
            grid_err=1
            goto 888
        endif
        
        ! panel side centers
        cx1=(x1+x2)/2
        cx2=(x2+x3)/2
        cx3=(x3+x4)/2
        cx4=(x4+x1)/2
        cy1=(y1+y2)/2
        cy2=(y2+y3)/2
        cy3=(y3+y4)/2
        cy4=(y4+y1)/2
        cz1=(z1+z2)/2
        cz2=(z2+z3)/2
        cz3=(z3+z4)/2
        cz4=(z4+z1)/2
        
        ! center of mass of quadrilateral panel
        cx=(cx1*d1 + cx2*d2 + cx3*d3 + cx4*d4)/dsum
        cy=(cy1*d1 + cy2*d2 + cy3*d3 + cy4*d4)/dsum
        cz=(cz1*d1 + cz2*d2 + cz3*d3 + cz4*d4)/dsum
        
        ! calculation of panel coordinates (in local CS)
        xl1=(x1-cx)*l1 + (y1-cy)*l2 + (z1-cz)*l3
        xl2=(x2-cx)*l1 + (y2-cy)*l2 + (z2-cz)*l3
        xl3=(x3-cx)*l1 + (y3-cy)*l2 + (z3-cz)*l3
        xl4=(x4-cx)*l1 + (y4-cy)*l2 + (z4-cz)*l3
        yl1=(x1-cx)*p1 + (y1-cy)*p2 + (z1-cz)*p3
        yl2=(x2-cx)*p1 + (y2-cy)*p2 + (z2-cz)*p3
        yl3=(x3-cx)*p1 + (y3-cy)*p2 + (z3-cz)*p3
        yl4=(x4-cx)*p1 + (y4-cy)*p2 + (z4-cz)*p3
    endif

    ! collocation points lowered by colldist
    colx=cx - colldist*n1
    coly=cy - colldist*n2
    colz=cz - colldist*n3

888 continue

    return
    end subroutine grid_panel_quad
    
! SUBROUTINE GRID_PANEL_TRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! subroutine for calculating grid data of triangular panels
    
    subroutine grid_panel_tri( interactive                                  ,&
                             & collcalc                                     ,&
                             & farfield                                     ,&
                             & error                                        ,&
                             & colldist                                     ,&
                             & x1,x2,x3                                     ,&
                             & y1,y2,y3                                     ,&
                             & z1,z2,z3                                     ,&
                             & FF                                           ,&
                             & S                                            ,&
                             & n1,n2,n3                                     ,&
                             & l1,l2,l3                                     ,&
                             & p1,p2,p3                                     ,&
                             & cx,cy,cz                                     ,&
                             & colx,coly,colz                               ,&
                             & xl1,xl2,xl3                                  ,&
                             & yl1,yl2,yl3                                  ,&
                             & d1,d2,d3                                     ,&
                             & grid_err                                     )

    implicit none
    
    ! INTENTS IN ===================================================================
    integer, intent(in)                    :: interactive                                    ,&
                                            & collcalc
    real(kind=kind_float),    intent(in)   :: x1,x2,x3                                       ,&
                                            & y1,y2,y3                                       ,&
                                            & z1,z2,z3                                       ,&
                                            & farfield                                       ,&
                                            & error                                          ,&
                                            & colldist
    ! INTENTS OUT ==================================================================
    integer, intent(out)                   :: grid_err
    real(kind=kind_float),    intent(out)  :: FF                                             ,&
                                            & S                                              ,&
                                            & n1,n2,n3                                       ,&
                                            & l1,l2,l3                                       ,&
                                            & p1,p2,p3                                       ,&
                                            & cx,cy,cz                                       ,&
                                            & colx,coly,colz                                 ,&
                                            & xl1,xl2,xl3                                    ,&
                                            & yl1,yl2,yl3                                    ,&
                                            & d1,d2,d3
    ! PRIVATE ======================================================================
    real(kind=kind_float)                  :: A_mod,B_mod,C_mod                              ,&
                                            & lx,ly,lz                                       ,&
                                            & cx1,cx2,cx3                                    ,&
                                            & cy1,cy2,cy3                                    ,&
                                            & cz1,cz2,cz3                                    ,&
                                            & l_mod                                          ,&
                                            & dsum
    real(kind=kind_float),    dimension(3) :: A,B,C,p
    ! ==============================================================================
    
    grid_err=0

    ! determining diagonals as vectors
    A=(/(x3-x1), (y3-y1), (z3-z1)/)
    B=(/(x3-x2), (y3-y2), (z3-z2)/)
    
    ! calculating their modules
    call func_vect_mod(A_mod,A)
    call func_vect_mod(B_mod,B)
    
    ! calculating far field distance as longer diagonal multiplied by farfield factor
    FF=farfield*max(A_mod,B_mod)
    
    ! calculating vector normal to diagonals (vector product)
    call func_vect_prod(C,A,B)
    
    ! calculating modulus of that vector
    call func_vect_mod(C_mod,C)
    
    ! panel surface is exactly half of that modulus
    S=C_mod/2._kind_float
    if (C_mod .lt. error) then
        call func_message( interactive                                      ,&
                         & "ERROR: Panel with zero surface detected")
        grid_err=1
        goto 777
    endif
    
    ! normalizing normal vector to unit vector
    n1=C(1)/C_mod
    n2=C(2)/C_mod
    n3=C(3)/C_mod
    
    ! calculating lengthwise unit vector (/l1,l2,l3/)
    lx=x3-(x1+x2)/2
    ly=y3-(y1+y2)/2
    lz=z3-(z1+z2)/2
    call func_vect_mod(l_mod,(/lx,ly,lz/))
    if (l_mod .lt. error) then
        call func_message( interactive                                      ,&
                         & "ERROR: Zero length lengthwise unit vector detected")
        grid_err=1
        goto 777
    endif
    l1=lx/l_mod
    l2=ly/l_mod
    l3=lz/l_mod
    
    ! calculating unit vector perpendicular to normal and lengthwise unit vectors (/p1,p2,p3/)
    call func_vect_prod(p,(/n1,n2,n3/),(/l1,l2,l3/))
    p1=p(1)
    p2=p(2)
    p3=p(3)
    
    if (collcalc .eq. 0) then
        ! approximative collocation/center points calculation
        
        ! center points
        cx=(x1 + x2 + x3)/3
        cy=(y1 + y2 + y3)/3
        cz=(z1 + z2 + z3)/3
        
        ! panel node coordinates (in local CS)
        xl1=(x1-cx)*l1 + (y1-cy)*l2 + (z1-cz)*l3
        xl2=(x2-cx)*l1 + (y2-cy)*l2 + (z2-cz)*l3
        xl3=(x3-cx)*l1 + (y3-cy)*l2 + (z3-cz)*l3
        yl1=(x1-cx)*p1 + (y1-cy)*p2 + (z1-cz)*p3
        yl2=(x2-cx)*p1 + (y2-cy)*p2 + (z2-cz)*p3
        yl3=(x3-cx)*p1 + (y3-cy)*p2 + (z3-cz)*p3
        
        ! panel side lengths
        d1=sqrt((xl2-xl1)**2 + (yl2-yl1)**2)
        d2=sqrt((xl3-xl2)**2 + (yl3-yl2)**2)
        d3=sqrt((xl1-xl3)**2 + (yl1-yl3)**2)
    elseif (collcalc .eq. 1) then
        ! accurate collocation/center points calculation
        
        ! panel side lengths and their sum
        call func_vect_mod(d1,(/(x2-x1),(y2-y1),(z2-z1)/))
        call func_vect_mod(d2,(/(x3-x2),(y3-y2),(z3-z2)/))
        call func_vect_mod(d3,(/(x1-x3),(y1-y3),(z1-z3)/))
        dsum=d1 + d2 + d3
        if (dsum .lt. error) then
            call func_message( interactive                                      ,&
                             & "ERROR: Zero sum of panel sides detected")
            grid_err=1
            goto 777
        endif
        
        ! panel side centers
        cx1=(x1+x2)/2
        cx2=(x2+x3)/2
        cx3=(x3+x1)/2
        cy1=(y1+y2)/2
        cy2=(y2+y3)/2
        cy3=(y3+y1)/2
        cz1=(z1+z2)/2
        cz2=(z2+z3)/2
        cz3=(z3+z1)/2
        
        ! center of mass of quadrilateral panel
        cx=(cx1*d1 + cx2*d2 + cx3*d3)/dsum
        cy=(cy1*d1 + cy2*d2 + cy3*d3)/dsum
        cz=(cz1*d1 + cz2*d2 + cz3*d3)/dsum
        
        ! calculation of panel coordinates (in local CS)
        xl1=(x1-cx)*l1 + (y1-cy)*l2 + (z1-cz)*l3
        xl2=(x2-cx)*l1 + (y2-cy)*l2 + (z2-cz)*l3
        xl3=(x3-cx)*l1 + (y3-cy)*l2 + (z3-cz)*l3
        yl1=(x1-cx)*p1 + (y1-cy)*p2 + (z1-cz)*p3
        yl2=(x2-cx)*p1 + (y2-cy)*p2 + (z2-cz)*p3
        yl3=(x3-cx)*p1 + (y3-cy)*p2 + (z3-cz)*p3
    endif

    ! collocation points lowered by colldist
    colx=cx - colldist*n1
    coly=cy - colldist*n2
    colz=cz - colldist*n3

777 continue

    return
    end subroutine grid_panel_tri
    
end module module_grid
