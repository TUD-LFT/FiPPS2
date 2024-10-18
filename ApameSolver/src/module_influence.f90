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

! file module_influence.f90

! This module calculates influence coefficient matrices a (and b) based on
! grid information

module module_influence

use module_kind_and_konst

! DECLARATIONS =================================================================
integer                                            :: influence_err
real(kind=kind_float), allocatable, dimension(:,:) :: a,a_no_wake,b
! ==============================================================================

contains

! SUBROUTINE INFLUENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subroutine influence( method                                            ,&
                        & symmetry                                          ,&
                        & panel_num                                         ,&
                        & panel_num_no_wake                                 ,&
                        & interactive                                       ,&
                        & panel_type                                        ,&
                        & elem1,elem2                                       ,&
                        & FF                                                ,&
                        & S                                                 ,&
                        & error                                             ,&
                        & colx,coly,colz                                    ,&
                        & cx,cy,cz                                          ,&
                        & l1,l2,l3                                          ,&
                        & p1,p2,p3                                          ,&
                        & n1,n2,n3                                          ,&
                        & x1,x2,x3,x4                                       ,&
                        & y1,y2,y3,y4                                       ,&
                        & d1,d2,d3,d4                                       ,&
                        & wake_only                                         )
    
    implicit none
    
    ! INTENTS IN ===================================================================
    integer,                       intent(in)                  :: method                     ,&
                                                                & symmetry                   ,&
                                                                & panel_num                  ,&
                                                                & panel_num_no_wake          ,&
                                                                & interactive
    real(kind=kind_float),                          intent(in) :: error
    integer, dimension(panel_num), intent(in)                  :: panel_type,elem1,elem2
    real(kind=kind_float),    dimension(panel_num), intent(in) :: FF                         ,&
                                                                & S                          ,&
                                                                & colx,coly,colz             ,&
                                                                & cx,cy,cz                   ,&
                                                                & l1,l2,l3                   ,&
                                                                & p1,p2,p3                   ,&
                                                                & n1,n2,n3                   ,&
                                                                & x1,x2,x3,x4                ,&
                                                                & y1,y2,y3,y4                ,&
                                                                & d1,d2,d3,d4
    logical,                                        intent(in) :: wake_only
    ! PRIVATE ======================================================================
    integer                                                    :: i,j                        ,&
                                                                & alloc_stat                 ,&
                                                                & nthreads                   ,&
                                                                & tid                        ,&
                                                                & omp_get_num_threads        ,&
                                                                & omp_get_thread_num         ,&
                                                                & chunksize                  ,&
                                                                & chunk
    real(kind=kind_float)                                      :: a_dummy, b_dummy
    parameter (chunksize=10)
    ! ==============================================================================

    ! getting number of threads from environment variable
!$  tid=omp_get_thread_num()
!$  if (tid .eq. 0) then
!$      nthreads=omp_get_num_threads()
!$  end if
    
    ! initialize "no errors" variable
    influence_err=0

    if (wake_only .eqv. .false.) then
        if (method .eq. 0) then
            ! constant source/doublet method
        
            ! allocate doublet influence matrix
            allocate ( a(panel_num_no_wake, panel_num_no_wake)         ,&
                    & a_no_wake(panel_num_no_wake, panel_num_no_wake)  ,&
                    & stat=alloc_stat                                  )
            if (alloc_stat .ne. 0) then
                call func_message( interactive                                          ,&
                    & "    ERROR: Not enough memory for allocating doublet influence coefficient matrix")
                influence_err=1
                goto 999
            endif
            
            ! allocate source influence matrix
            allocate ( b(panel_num_no_wake, panel_num_no_wake)  ,&
                     & stat=alloc_stat                          )
            if (alloc_stat .ne. 0) then
                call func_message( interactive                                          ,&
                    & "    ERROR: Not enough memory for allocating source influence coefficient matrix,&
                    & consider using doublet only singularity method")
                influence_err=1
                goto 999
            endif
            
            chunk=chunksize
            !$OMP PARALLEL SHARED(panel_type,a_no_wake,b,FF,S,colx,coly,colz,cx,cy,cz,l1,l2,l3,&
            !$OMP &p1,p2,p3,n1,n2,n3,x1,x2,x3,x4,y1,y2,y3,y4,d1,d2,d3,d4,nthreads,chunk) PRIVATE(i,tid,a_dummy,b_dummy)
            !$OMP DO SCHEDULE(static,chunk)
            
            do i=1,panel_num_no_wake        ! influenced panel
                do j=1,panel_num_no_wake    ! influencing panel
                    ! calculate influence coefficient
                    if (panel_type(j) .eq. 1 .or. panel_type(j) .eq. 20) then
                        ! quadrilateral panel
                        call influence_src_dub_panel_quad( a_no_wake(i,j),b(i,j)    ,&
                                                         & FF(j)                    ,&
                                                         & S(j)                     ,&
                                                         & error                    ,&
                                                         & colx(i),coly(i),colz(i)  ,&
                                                         & cx(j),cy(j),cz(j)        ,&
                                                         & l1(j),l2(j),l3(j)        ,&
                                                         & p1(j),p2(j),p3(j)        ,&
                                                         & n1(j),n2(j),n3(j)        ,&
                                                         & x1(j),x2(j),x3(j),x4(j)  ,&
                                                         & y1(j),y2(j),y3(j),y4(j)  ,&
                                                         & d1(j),d2(j),d3(j),d4(j)  )
                    elseif (panel_type(j) .eq. 2 .or. panel_type(j) .eq. 21) then
                        ! triangular panel
                        call influence_src_dub_panel_tri( a_no_wake(i,j),b(i,j)     ,&
                                                        & FF(j)                     ,&
                                                        & S(j)                      ,&
                                                        & error                     ,&
                                                        & colx(i),coly(i),colz(i)   ,&
                                                        & cx(j),cy(j),cz(j)         ,&
                                                        & l1(j),l2(j),l3(j)         ,&
                                                        & p1(j),p2(j),p3(j)         ,&
                                                        & n1(j),n2(j),n3(j)         ,&
                                                        & x1(j),x2(j),x3(j)         ,&
                                                        & y1(j),y2(j),y3(j)         ,&
                                                        & d1(j),d2(j),d3(j)         )
                    endif
                    
                    ! self influence
                    if (i .eq. j) then
                        a_no_wake(i,j)=-0.5_kind_float
                    endif
                    ! add symmetry conditions (symmetry about x-z-plane)
                    if (symmetry .eq. 1) then
                        ! calculate influence coefficient
                        if (panel_type(j) .eq. 1 .or. panel_type(j) .eq. 20) then
                            ! quadrilateral panel
                            call influence_src_dub_panel_quad( a_dummy,b_dummy         ,&
                                                            & FF(j)                    ,&
                                                            & S(j)                     ,&
                                                            & error                    ,&
                                                            & colx(i),-coly(i),colz(i) ,&
                                                            & cx(j),cy(j),cz(j)        ,&
                                                            & l1(j),l2(j),l3(j)        ,&
                                                            & p1(j),p2(j),p3(j)        ,&
                                                            & n1(j),n2(j),n3(j)        ,&
                                                            & x1(j),x2(j),x3(j),x4(j)  ,&
                                                            & y1(j),y2(j),y3(j),y4(j)  ,&
                                                            & d1(j),d2(j),d3(j),d4(j)  )
                        elseif (panel_type(j) .eq. 2 .or. panel_type(j) .eq. 21) then
                            ! triangular panel
                            call influence_src_dub_panel_tri( a_dummy,b_dummy           ,&
                                                            & FF(j)                     ,&
                                                            & S(j)                      ,&
                                                            & error                     ,&
                                                            & colx(i),-coly(i),colz(i)  ,&
                                                            & cx(j),cy(j),cz(j)         ,&
                                                            & l1(j),l2(j),l3(j)         ,&
                                                            & p1(j),p2(j),p3(j)         ,&
                                                            & n1(j),n2(j),n3(j)         ,&
                                                            & x1(j),x2(j),x3(j)         ,&
                                                            & y1(j),y2(j),y3(j)         ,&
                                                            & d1(j),d2(j),d3(j)         )
                        endif
                        a_no_wake(i,j) = a_no_wake(i,j) + a_dummy
                        b(i,j)         = b(i,j)         + b_dummy
                    end if
                enddo
            enddo
            
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

        elseif (method .eq. 1) then
            ! constant doublet method
        
            ! allocate only a (doublet) matrix
            allocate ( a(panel_num_no_wake, panel_num_no_wake)          ,&
                     & a_no_wake(panel_num_no_wake, panel_num_no_wake)  ,&
                     & stat=alloc_stat                                  )
            if (alloc_stat .ne. 0) then
                call func_message( interactive                                          ,&
                    & "    ERROR: Not enough memory for allocating influence coefficient matrix")
                influence_err=1
                goto 999
            endif
            
            chunk=chunksize
            
            !$OMP PARALLEL SHARED(panel_type,a_no_wake,FF,S,colx,coly,colz,cx,cy,cz,l1,l2,l3,p1,p2,p3,&
            !$OMP &n1,n2,n3,x1,x2,x3,x4,y1,y2,y3,y4,d1,d2,d3,d4,nthreads,chunk) PRIVATE(i,tid,a_dummy)
            !$OMP DO SCHEDULE(static,chunk)
            
            do i=1,panel_num_no_wake        ! influenced panel
                do j=1,panel_num_no_wake    ! influencing panel
                    ! calculate influence coefficient
                    if (i .eq. j) then
                        ! self influence
                        a_no_wake(i,j)=-0.5_kind_float
                    elseif (panel_type(j) .eq. 1 .or. panel_type(j) .eq. 20) then
                        ! quadrilateral panel
                        call influence_dub_panel_quad( a_no_wake(i,j)           ,&
                                                     & FF(j)                    ,&
                                                     & S(j)                     ,&
                                                     & error                    ,&
                                                     & colx(i),coly(i),colz(i)  ,&
                                                     & cx(j),cy(j),cz(j)        ,&
                                                     & l1(j),l2(j),l3(j)        ,&
                                                     & p1(j),p2(j),p3(j)        ,&
                                                     & n1(j),n2(j),n3(j)        ,&
                                                     & x1(j),x2(j),x3(j),x4(j)  ,&
                                                     & y1(j),y2(j),y3(j),y4(j)  ,&
                                                     & d1(j),d2(j),d3(j),d4(j)  )
                    elseif (panel_type(j) .eq. 2 .or. panel_type(j) .eq. 21) then
                        ! triangular panel
                        call influence_dub_panel_tri( a_no_wake(i,j)            ,&
                                                    & FF(j)                     ,&
                                                    & S(j)                      ,&
                                                    & error                     ,&
                                                    & colx(i),coly(i),colz(i)   ,&
                                                    & cx(j),cy(j),cz(j)         ,&
                                                    & l1(j),l2(j),l3(j)         ,&
                                                    & p1(j),p2(j),p3(j)         ,&
                                                    & n1(j),n2(j),n3(j)         ,&
                                                    & x1(j),x2(j),x3(j)         ,&
                                                    & y1(j),y2(j),y3(j)         ,&
                                                    & d1(j),d2(j),d3(j)         )
                    endif
                    ! add symmetry conditions to (symmetry about x-z-plane)
                    if (symmetry .eq. 1) then
                        if (panel_type(j) .eq. 1 .or. panel_type(j) .eq. 20) then
                            ! quadrilateral panel
                            call influence_dub_panel_quad( a_dummy                 ,&
                                                        & FF(j)                    ,&
                                                        & S(j)                     ,&
                                                        & error                    ,&
                                                        & colx(i),-coly(i),colz(i) ,&
                                                        & cx(j),cy(j),cz(j)        ,&
                                                        & l1(j),l2(j),l3(j)        ,&
                                                        & p1(j),p2(j),p3(j)        ,&
                                                        & n1(j),n2(j),n3(j)        ,&
                                                        & x1(j),x2(j),x3(j),x4(j)  ,&
                                                        & y1(j),y2(j),y3(j),y4(j)  ,&
                                                        & d1(j),d2(j),d3(j),d4(j)  )
                        elseif (panel_type(j) .eq. 2 .or. panel_type(j) .eq. 21) then
                            ! triangular panel
                            call influence_dub_panel_tri( a_dummy                   ,&
                                                        & FF(j)                     ,&
                                                        & S(j)                      ,&
                                                        & error                     ,&
                                                        & colx(i),-coly(i),colz(i)  ,&
                                                        & cx(j),cy(j),cz(j)         ,&
                                                        & l1(j),l2(j),l3(j)         ,&
                                                        & p1(j),p2(j),p3(j)         ,&
                                                        & n1(j),n2(j),n3(j)         ,&
                                                        & x1(j),x2(j),x3(j)         ,&
                                                        & y1(j),y2(j),y3(j)         ,&
                                                        & d1(j),d2(j),d3(j)         )
                        end if
                        a_no_wake(i,j) = a_no_wake(i,j) + a_dummy
                    end if
                enddo
            enddo
            
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            
        endif
    end if

    a(:,:) = a_no_wake(:,:)

    ! now add influence of wake panels
    do i=1,panel_num_no_wake                ! influenced panel
        do j=panel_num_no_wake+1,panel_num  ! influencing panel
            ! calculate influence coefficient
            if (panel_type(j) .eq. 10) then
                call influence_dub_panel_quad( a_dummy                  ,&
                                             & FF(j)                    ,&
                                             & S(j)                     ,&
                                             & error                    ,&
                                             & colx(i),coly(i),colz(i)  ,&
                                             & cx(j),cy(j),cz(j)        ,&
                                             & l1(j),l2(j),l3(j)        ,&
                                             & p1(j),p2(j),p3(j)        ,&
                                             & n1(j),n2(j),n3(j)        ,&
                                             & x1(j),x2(j),x3(j),x4(j)  ,&
                                             & y1(j),y2(j),y3(j),y4(j)  ,&
                                             & d1(j),d2(j),d3(j),d4(j)  )
                a(i,(elem1(j))) = a(i,(elem1(j))) + a_dummy
                a(i,(elem2(j))) = a(i,(elem2(j))) - a_dummy
            elseif (panel_type(j) .eq. 11) then
                ! triangular wake panel
                call influence_dub_panel_tri( a_dummy                   ,&
                                            & FF(j)                     ,&
                                            & S(j)                      ,&
                                            & error                     ,&
                                            & colx(i),coly(i),colz(i)   ,&
                                            & cx(j),cy(j),cz(j)         ,&
                                            & l1(j),l2(j),l3(j)         ,&
                                            & p1(j),p2(j),p3(j)         ,&
                                            & n1(j),n2(j),n3(j)         ,&
                                            & x1(j),x2(j),x3(j)         ,&
                                            & y1(j),y2(j),y3(j)         ,&
                                            & d1(j),d2(j),d3(j)         )
                a(i,(elem1(j))) = a(i,(elem1(j))) + a_dummy
                a(i,(elem2(j))) = a(i,(elem2(j))) - a_dummy
            endif
            if (symmetry .eq. 1) then
                ! calculate influence coefficient
                if (panel_type(j) .eq. 10) then
                    call influence_dub_panel_quad( a_dummy                  ,&
                                                & FF(j)                    ,&
                                                & S(j)                     ,&
                                                & error                    ,&
                                                & colx(i),-coly(i),colz(i) ,&
                                                & cx(j),cy(j),cz(j)        ,&
                                                & l1(j),l2(j),l3(j)        ,&
                                                & p1(j),p2(j),p3(j)        ,&
                                                & n1(j),n2(j),n3(j)        ,&
                                                & x1(j),x2(j),x3(j),x4(j)  ,&
                                                & y1(j),y2(j),y3(j),y4(j)  ,&
                                                & d1(j),d2(j),d3(j),d4(j)  )
                    a(i,(elem1(j))) = a(i,(elem1(j))) + a_dummy
                    a(i,(elem2(j))) = a(i,(elem2(j))) - a_dummy
                elseif (panel_type(j) .eq. 11) then
                    ! triangular wake panel
                    call influence_dub_panel_tri( a_dummy                   ,&
                                                & FF(j)                     ,&
                                                & S(j)                      ,&
                                                & error                     ,&
                                                & colx(i),-coly(i),colz(i)  ,&
                                                & cx(j),cy(j),cz(j)         ,&
                                                & l1(j),l2(j),l3(j)         ,&
                                                & p1(j),p2(j),p3(j)         ,&
                                                & n1(j),n2(j),n3(j)         ,&
                                                & x1(j),x2(j),x3(j)         ,&
                                                & y1(j),y2(j),y3(j)         ,&
                                                & d1(j),d2(j),d3(j)         )
                    a(i,(elem1(j))) = a(i,(elem1(j))) + a_dummy
                    a(i,(elem2(j))) = a(i,(elem2(j))) - a_dummy
                endif
            end if
        enddo
    enddo
    
999 continue
    
    return
    end subroutine influence
    
! SUBROUTINE INFLUENCE_SRC_DUB_PANEL_QUAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine calculates influence coefficient of quadrilateral panel
    ! for constant source/dipole distribution
    
    subroutine influence_src_dub_panel_quad( a,b                            ,&
                                           & FF                             ,&
                                           & S                              ,&
                                           & error                          ,&
                                           & colx,coly,colz                 ,&
                                           & cx,cy,cz                       ,&
                                           & l1,l2,l3                       ,&
                                           & p1,p2,p3                       ,&
                                           & n1,n2,n3                       ,&
                                           & x1,x2,x3,x4                    ,&
                                           & y1,y2,y3,y4                    ,&
                                           & d1,d2,d3,d4                    )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float), intent(in)  :: FF                                                 ,&
                                        & S                                                  ,&
                                        & error                                              ,&
                                        & colx,coly,colz                                     ,&
                                        & cx,cy,cz                                           ,&
                                        & l1,l2,l3                                           ,&
                                        & p1,p2,p3                                           ,&
                                        & n1,n2,n3                                           ,&
                                        & x1,x2,x3,x4                                        ,&
                                        & y1,y2,y3,y4                                        ,&
                                        & d1,d2,d3,d4
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), intent(out) :: a,b
    ! PRIVATE ======================================================================
    real(kind=kind_float)              :: fourpi                                             ,&
                                        & FF_sqr                                             ,&
                                        & dist_x,dist_y,dist_z                               ,&
                                        & cpx,cpy,cpz                                        ,&
                                        & rad                                                ,&
                                        & cpx1,cpx2,cpx3,cpx4                                ,&
                                        & cpy1,cpy2,cpy3,cpy4                                ,&
                                        & e1,e2,e3,e4                                        ,&
                                        & r1,r2,r3,r4                                        ,&
                                        & h1,h2,h3,h4                                        ,&
                                        & x21,x32,x43,x14                                    ,&
                                        & y21,y32,y43,y14                                    ,&
                                        & F                                                  ,&
                                        & G                                                  ,&
                                        & a1,a2,a3,a4                                        ,&
                                        & b1,b2,b3,b4
    parameter (fourpi=0.079577471545948_kind_float)
    ! ==============================================================================

    FF_sqr=FF**2
    dist_x=colx-cx
    dist_y=coly-cy
    dist_z=colz-cz
    cpx=dist_x*l1 + dist_y*l2 + dist_z*l3
    cpy=dist_x*p1 + dist_y*p2 + dist_z*p3
    cpz=dist_x*n1 + dist_y*n2 + dist_z*n3
    rad=cpx**2 + cpy**2 + cpz**2
    ! if distance of panel from influenced point is greater
    ! then product of longer diagonal and "far field" coefficient
    if (rad .gt. FF_sqr) then
        a=fourpi*S*cpz*rad**(-1.5_kind_float)
        b=fourpi*S/sqrt(rad)
    else
        cpx1=cpx - x1
        cpx2=cpx - x2
        cpx3=cpx - x3
        cpx4=cpx - x4
        cpy1=cpy - y1
        cpy2=cpy - y2
        cpy3=cpy - y3
        cpy4=cpy - y4
        e1=cpx1**2+cpz**2
        e2=cpx2**2+cpz**2
        e3=cpx3**2+cpz**2
        e4=cpx4**2+cpz**2
        r1=sqrt(e1 + cpy1**2)
        r2=sqrt(e2 + cpy2**2)
        r3=sqrt(e3 + cpy3**2)
        r4=sqrt(e4 + cpy4**2)
        x21=x2-x1
        x32=x3-x2
        x43=x4-x3
        x14=x1-x4
        y21=y2-y1
        y32=y3-y2
        y43=y4-y3
        y14=y1-y4
        if (abs(cpz) .lt. error) then
            if (d1 .lt. error) then
                b1=0
            else
                b1=(cpx1*y21-cpy1*x21)/d1*log((r1+r2+d1)/(r1+r2-d1))
            endif
            if (d2 .lt. error) then
                b2=0
            else
                b2=(cpx2*y32-cpy2*x32)/d2*log((r2+r3+d2)/(r2+r3-d2))
            endif
            if (d3 .lt. error) then
                b3=0
            else
                b3=(cpx3*y43-cpy3*x43)/d3*log((r3+r4+d3)/(r3+r4-d3))
            endif
            if (d4 .lt. error) then
                b4=0
            else
                b4=(cpx4*y14-cpy4*x14)/d4*log((r4+r1+d4)/(r4+r1-d4))
            endif
            a=0
            b=-(b1+b2+b3+b4)*fourpi
        else
            h1=cpx1*cpy1
            h2=cpx2*cpy2
            h3=cpx3*cpy3
            h4=cpx4*cpy4
            if (d1 .lt. error) then
                a1=0
                b1=0
            else
                F=y21*e1 - x21*h1
                G=y21*e2 - x21*h2
                a1=atan2(cpz*x21*(F*r2-G*r1), cpz**2*x21**2*r1*r2+F*G)
                b1=(cpx1*y21-cpy1*x21)/d1*log((r1+r2+d1)/(r1+r2-d1))
            endif
            if (d2 .lt. error) then
                a2=0
                b2=0
            else
                F=y32*e2 - x32*h2
                G=y32*e3 - x32*h3
                a2=atan2(cpz*x32*(F*r3-G*r2), cpz**2*x32**2*r2*r3+F*G)
                b2=(cpx2*y32-cpy2*x32)/d2*log((r2+r3+d2)/(r2+r3-d2))
            endif
            if (d3 .lt. error) then
                a3=0
                b3=0
            else
                F=y43*e3 - x43*h3
                G=y43*e4 - x43*h4
                a3=atan2(cpz*x43*(F*r4-G*r3), cpz**2*x43**2*r3*r4+F*G)
                b3=(cpx3*y43-cpy3*x43)/d3*log((r3+r4+d3)/(r3+r4-d3))
            endif
            if (d4 .lt. error) then
                a4=0
                b4=0
            else
                F=y14*e4 - x14*h4
                G=y14*e1 - x14*h1
                a4=atan2(cpz*x14*(F*r1-G*r4), cpz**2*x14**2*r4*r1+F*G)
                b4=(cpx4*y14-cpy4*x14)/d4*log((r4+r1+d4)/(r4+r1-d4))
            endif
            a=-(a1+a2+a3+a4)*fourpi
            b=-(b1+b2+b3+b4)*fourpi-cpz*a
        endif
    endif

    return
    end subroutine influence_src_dub_panel_quad
    
! SUBROUTINE INFLUENCE_SRC_DUB_PANEL_TRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine calculates influence coefficient of triangular panel
    ! for constant source/dipole distribution

    subroutine influence_src_dub_panel_tri( a,b                             ,&
                                          & FF                              ,&
                                          & S                               ,&
                                          & error                           ,&
                                          & colx,coly,colz                  ,&
                                          & cx,cy,cz                        ,&
                                          & l1,l2,l3                        ,&
                                          & p1,p2,p3                        ,&
                                          & n1,n2,n3                        ,&
                                          & x1,x2,x3                        ,&
                                          & y1,y2,y3                        ,&
                                          & d1,d2,d3                        )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float), intent(in)  :: FF                                                 ,&
                                        & S                                                  ,&
                                        & error                                              ,&
                                        & colx,coly,colz                                     ,&
                                        & cx,cy,cz                                           ,&
                                        & l1,l2,l3                                           ,&
                                        & p1,p2,p3                                           ,&
                                        & n1,n2,n3                                           ,&
                                        & x1,x2,x3                                           ,&
                                        & y1,y2,y3                                           ,&
                                        & d1,d2,d3
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), intent(out) :: a,b
    ! PRIVATE ======================================================================
    real(kind=kind_float)              :: fourpi                                             ,&
                                        & FF_sqr                                             ,&
                                        & dist_x,dist_y,dist_z                               ,&
                                        & cpx,cpy,cpz                                        ,&
                                        & rad                                                ,&
                                        & cpx1,cpx2,cpx3                                     ,&
                                        & cpy1,cpy2,cpy3                                     ,&
                                        & e1,e2,e3                                           ,&
                                        & r1,r2,r3                                           ,&
                                        & h1,h2,h3                                           ,&
                                        & x21,x32,x13                                        ,&
                                        & y21,y32,y13                                        ,&
                                        & F,G                                                ,&
                                        & a1,a2,a3                                           ,&
                                        & b1,b2,b3
    parameter (fourpi=0.079577471545948_kind_float)
    ! ==============================================================================
    
    FF_sqr=FF**2
    dist_x=colx-cx
    dist_y=coly-cy
    dist_z=colz-cz
    cpx=dist_x*l1 + dist_y*l2 + dist_z*l3
    cpy=dist_x*p1 + dist_y*p2 + dist_z*p3
    cpz=dist_x*n1 + dist_y*n2 + dist_z*n3
    rad=cpx**2 + cpy**2 + cpz**2
    ! if distance of panel from influenced point is greater
    ! then product of longer diagonal and "far field" coefficient
    if (rad .gt. FF_sqr) then
        a=fourpi*S*cpz*rad**(-1.5_kind_float)
        b=fourpi*S/sqrt(rad)
    else
        cpx1=cpx - x1
        cpx2=cpx - x2
        cpx3=cpx - x3
        cpy1=cpy - y1
        cpy2=cpy - y2
        cpy3=cpy - y3
        e1=cpx1**2+cpz**2
        e2=cpx2**2+cpz**2
        e3=cpx3**2+cpz**2
        r1=sqrt(e1 + cpy1**2)
        r2=sqrt(e2 + cpy2**2)
        r3=sqrt(e3 + cpy3**2)
        x21=x2-x1
        x32=x3-x2
        x13=x1-x3
        y21=y2-y1
        y32=y3-y2
        y13=y1-y3
        if (abs(cpz) .lt. error) then
            if (d1 .lt. error) then
                b1=0
            else
                b1=(cpx1*y21-cpy1*x21)/d1*log((r1+r2+d1)/(r1+r2-d1))
            endif
            if (d2 .lt. error) then
                b2=0
            else
                b2=(cpx2*y32-cpy2*x32)/d2*log((r2+r3+d2)/(r2+r3-d2))
            endif
            if (d3 .lt. error) then
                b3=0
            else
                b3=(cpx3*y13-cpy3*x13)/d3*log((r3+r1+d3)/(r3+r1-d3))
            endif
            a=0
            b=-(b1+b2+b3)*fourpi
        else
            h1=cpx1*cpy1
            h2=cpx2*cpy2
            h3=cpx3*cpy3
            if (d1 .lt. error) then
                a1=0
                b1=0
            else
                F=y21*e1 - x21*h1
                G=y21*e2 - x21*h2
                a1=atan2(cpz*x21*(F*r2-G*r1), cpz**2*x21**2*r1*r2+F*G)
                b1=(cpx1*y21-cpy1*x21)/d1*log((r1+r2+d1)/(r1+r2-d1))
            endif
            if (d2 .lt. error) then
                a2=0
                b2=0
            else
                F=y32*e2 - x32*h2
                G=y32*e3 - x32*h3
                a2=atan2(cpz*x32*(F*r3-G*r2), cpz**2*x32**2*r2*r3+F*G)
                b2=(cpx2*y32-cpy2*x32)/d2*log((r2+r3+d2)/(r2+r3-d2))
            endif
            if (d3 .lt. error) then
                a3=0
                b3=0
            else
                F=y13*e3 - x13*h3
                G=y13*e1 - x13*h1
                a3=atan2(cpz*x13*(F*r1-G*r3), cpz**2*x13**2*r3*r1+F*G)
                b3=(cpx3*y13-cpy3*x13)/d3*log((r3+r1+d3)/(r3+r1-d3))
            endif
            a=-(a1+a2+a3)*fourpi
            b=-(b1+b2+b3)*fourpi-cpz*a
        endif
    endif

    return
    end subroutine influence_src_dub_panel_tri
    
! SUBROUTINE INFLUENCE_DUB_PANEL_QUAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine calculates influence coefficient of quadrilateral panel
    ! for constant dipole distribution

    subroutine influence_dub_panel_quad( a                                  ,&
                                       & FF                                 ,&
                                       & S                                  ,&
                                       & error                              ,&
                                       & colx,coly,colz                     ,&
                                       & cx,cy,cz                           ,&
                                       & l1,l2,l3                           ,&
                                       & p1,p2,p3                           ,&
                                       & n1,n2,n3                           ,&
                                       & x1,x2,x3,x4                        ,&
                                       & y1,y2,y3,y4                        ,&
                                       & d1,d2,d3,d4                        )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float), intent(in)  :: FF                                                 ,&
                                        & S                                                  ,&
                                        & error                                              ,&
                                        & colx,coly,colz                                     ,&
                                        & cx,cy,cz                                           ,&
                                        & l1,l2,l3                                           ,&
                                        & p1,p2,p3                                           ,&
                                        & n1,n2,n3                                           ,&
                                        & x1,x2,x3,x4                                        ,&
                                        & y1,y2,y3,y4                                        ,&
                                        & d1,d2,d3,d4
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), intent(out) :: a
    ! PRIVATE ======================================================================
    real(kind=kind_float)              :: fourpi                                             ,&
                                        & FF_sqr                                             ,&
                                        & dist_x,dist_y,dist_z                               ,&
                                        & cpx,cpy,cpz                                        ,&
                                        & rad                                                ,&
                                        & cpx1,cpx2,cpx3,cpx4                                ,&
                                        & cpy1,cpy2,cpy3,cpy4                                ,&
                                        & e1,e2,e3,e4                                        ,&
                                        & r1,r2,r3,r4                                        ,&
                                        & h1,h2,h3,h4                                        ,&
                                        & x21,x32,x43,x14                                    ,&
                                        & y21,y32,y43,y14                                    ,&
                                        & F,G                                                ,&
                                        & a1,a2,a3,a4
    parameter (fourpi=0.079577471545948_kind_float)
    ! ==============================================================================
    
    FF_sqr=FF**2
    dist_x=colx-cx
    dist_y=coly-cy
    dist_z=colz-cz
    cpx=dist_x*l1 + dist_y*l2 + dist_z*l3
    cpy=dist_x*p1 + dist_y*p2 + dist_z*p3
    cpz=dist_x*n1 + dist_y*n2 + dist_z*n3
    rad=cpx**2 + cpy**2 + cpz**2
    ! if distance of panel from influenced point is greater
    ! then product of longer diagonal and "far field" coefficient
    if (rad .gt. FF_sqr) then
        a=fourpi*S*cpz*rad**(-1.5_kind_float)
    else
        cpx1=cpx - x1
        cpx2=cpx - x2
        cpx3=cpx - x3
        cpx4=cpx - x4
        cpy1=cpy - y1
        cpy2=cpy - y2
        cpy3=cpy - y3
        cpy4=cpy - y4
        e1=cpx1**2+cpz**2
        e2=cpx2**2+cpz**2
        e3=cpx3**2+cpz**2
        e4=cpx4**2+cpz**2
        r1=sqrt(e1 + cpy1**2)
        r2=sqrt(e2 + cpy2**2)
        r3=sqrt(e3 + cpy3**2)
        r4=sqrt(e4 + cpy4**2)
        x21=x2-x1
        x32=x3-x2
        x43=x4-x3
        x14=x1-x4
        y21=y2-y1
        y32=y3-y2
        y43=y4-y3
        y14=y1-y4
        if (abs(cpz) .lt. error) then
            a=0
        else
            h1=cpx1*cpy1
            h2=cpx2*cpy2
            h3=cpx3*cpy3
            h4=cpx4*cpy4
            if (d1 .lt. error) then
                a1=0
            else
                F=y21*e1 - x21*h1
                G=y21*e2 - x21*h2
                a1=atan2(cpz*x21*(F*r2-G*r1), cpz**2*x21**2*r1*r2+F*G)
            endif
            if (d2 .lt. error) then
                a2=0
            else
                F=y32*e2 - x32*h2
                G=y32*e3 - x32*h3
                a2=atan2(cpz*x32*(F*r3-G*r2), cpz**2*x32**2*r2*r3+F*G)
            endif
            if (d3 .lt. error) then
                a3=0
            else
                F=y43*e3 - x43*h3
                G=y43*e4 - x43*h4
                a3=atan2(cpz*x43*(F*r4-G*r3), cpz**2*x43**2*r3*r4+F*G)
            endif
            if (d4 .lt. error) then
                a4=0
            else
                F=y14*e4 - x14*h4
                G=y14*e1 - x14*h1
                a4=atan2(cpz*x14*(F*r1-G*r4), cpz**2*x14**2*r4*r1+F*G)
            endif
            a=-(a1+a2+a3+a4)*fourpi
        endif
    endif

    return
    end subroutine influence_dub_panel_quad
    
! SUBROUTINE INFLUENCE_DUB_PANEL_TRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine calculates influence coefficient of triangular panel
    ! for constant dipole distribution

    subroutine influence_dub_panel_tri( a                                   ,&
                                      & FF                                  ,&
                                      & S                                   ,&
                                      & error                               ,&
                                      & colx,coly,colz                      ,&
                                      & cx,cy,cz                            ,&
                                      & l1,l2,l3                            ,&
                                      & p1,p2,p3                            ,&
                                      & n1,n2,n3                            ,&
                                      & x1,x2,x3                            ,&
                                      & y1,y2,y3                            ,&
                                      & d1,d2,d3                            )

    implicit none

    ! INTENTS IN ===================================================================
    real(kind=kind_float), intent(in)  :: FF                                                 ,&
                                        & S                                                  ,&
                                        & error                                              ,&
                                        & colx,coly,colz                                     ,&
                                        & cx,cy,cz                                           ,&
                                        & l1,l2,l3                                           ,&
                                        & p1,p2,p3                                           ,&
                                        & n1,n2,n3                                           ,&
                                        & x1,x2,x3                                           ,&
                                        & y1,y2,y3                                           ,&
                                        & d1,d2,d3
    ! INTENTS OUT ==================================================================
    real(kind=kind_float), intent(out) :: a
    ! PRIVATE ======================================================================
    real(kind=kind_float)              :: fourpi                                             ,&
                                        & FF_sqr                                             ,&
                                        & dist_x,dist_y,dist_z                               ,&
                                        & cpx,cpy,cpz                                        ,&
                                        & rad                                                ,&
                                        & cpx1,cpx2,cpx3                                     ,&
                                        & cpy1,cpy2,cpy3                                     ,&
                                        & e1,e2,e3                                           ,&
                                        & r1,r2,r3                                           ,&
                                        & h1,h2,h3                                           ,&
                                        & x21,x32,x13                                        ,&
                                        & y21,y32,y13                                        ,&
                                        & F,G                                                ,&
                                        & a1,a2,a3
    parameter (fourpi=0.079577471545948_kind_float)
    ! ==============================================================================
    
    FF_sqr=FF**2
    dist_x=colx-cx
    dist_y=coly-cy
    dist_z=colz-cz
    cpx=dist_x*l1 + dist_y*l2 + dist_z*l3
    cpy=dist_x*p1 + dist_y*p2 + dist_z*p3
    cpz=dist_x*n1 + dist_y*n2 + dist_z*n3
    rad=cpx**2 + cpy**2 + cpz**2
    ! if distance of panel from influenced point is greater
    ! then product of longer diagonal and "far field" coefficient
    if (rad .gt. FF_sqr) then
        a=fourpi*S*cpz*rad**(-1.5_kind_float)
    else
        cpx1=cpx - x1
        cpx2=cpx - x2
        cpx3=cpx - x3
        cpy1=cpy - y1
        cpy2=cpy - y2
        cpy3=cpy - y3
        e1=cpx1**2+cpz**2
        e2=cpx2**2+cpz**2
        e3=cpx3**2+cpz**2
        r1=sqrt(e1 + cpy1**2)
        r2=sqrt(e2 + cpy2**2)
        r3=sqrt(e3 + cpy3**2)
        x21=x2-x1
        x32=x3-x2
        x13=x1-x3
        y21=y2-y1
        y32=y3-y2
        y13=y1-y3
        if (abs(cpz) .lt. error) then
            a=0
        else
            h1=cpx1*cpy1
            h2=cpx2*cpy2
            h3=cpx3*cpy3
            if (d1 .lt. error) then
                a1=0
            else
                F=y21*e1 - x21*h1
                G=y21*e2 - x21*h2
                a1=atan2(cpz*x21*(F*r2-G*r1), cpz**2*x21**2*r1*r2+F*G)
            endif
            if (d2 .lt. error) then
                a2=0
            else
                F=y32*e2 - x32*h2
                G=y32*e3 - x32*h3
                a2=atan2(cpz*x32*(F*r3-G*r2), cpz**2*x32**2*r2*r3+F*G)
            endif
            if (d3 .lt. error) then
                a3=0
            else
                F=y13*e3 - x13*h3
                G=y13*e1 - x13*h1
                a3=atan2(cpz*x13*(F*r1-G*r3), cpz**2*x13**2*r3*r1+F*G)
            endif
            a=-(a1+a2+a3)*fourpi
        endif
    endif

    return
    end subroutine influence_dub_panel_tri

end module module_influence
