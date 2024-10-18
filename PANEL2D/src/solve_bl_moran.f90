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
subroutine solve_bl_moran(verbose,Re,nodes,nblp,ue,dxsep,transloc,n_te,bl_values,phi_te,error)

! Load modules
use konst_var
use functions
use bl_functions

implicit none

! Declaration of variables
intent(in)                                        :: verbose,Re,nodes,nblp
intent(in)                                        :: ue
intent(out)                                       :: transloc,dxsep,bl_values,error,n_te,phi_te
integer                                           :: nblp,n_te
double precision                                  :: Re,dxsep,phi_te
double precision, dimension(nblp,2)               :: nodes
double precision, dimension(nblp)                 :: ue
type(bl_point), dimension(nblp)                   :: bl_values
logical                                           :: verbose

double precision                                  :: dth2ue6
double precision                                  :: lambda,L
double precision                                  :: rex,ret,retmax
double precision                                  :: transloc,tsep
double precision                                  :: rtheta
double precision                                  :: scrit,H_guess,H_delta
double precision, dimension(nblp)                 :: s
double precision, dimension(2)                    :: y, temp

double precision, dimension(nblp)                 :: dueds

integer                                           :: itrans,endofsurf,transorlamsep,turbseparation
integer                                           :: ii,jj,isep

logical                                           :: error

error=.false.

tsep           = 0.d0
turbseparation = 0

bl_values(:)%H      = 0.d0
bl_values(:)%theta  = 0.d0
bl_values(:)%deltas = 0.d0
bl_values(:)%cf     = 0.d0
bl_values(:)%nfac   = 0.d0

! Running coordinate along airfoil
call running_coordinate(nblp,nodes,s)

! Velocity gradients
call velocity_gradient(nblp,ue,s,dueds)

! Laminar boundary layer
! According to THWAITES
transorlamsep = 0
endofsurf     = 0
scrit         = 0.d0

bl_values(1)%theta = sqrt(0.075d0/Re/dueds(1))

ii = 1
loop1: do while((transorlamsep .eq. 0) .and. (endofsurf .eq. 0))

    if (ii .eq. 2) then
        bl_values(2)%theta = bl_values(1)%theta
    else if (ii .gt. 2) then
        dth2ue6   = 0.225d0*(ue(ii)**5 + ue(ii-1)**5)*(s(ii) - s(ii-1))/Re
        bl_values(ii)%theta = sqrt(((bl_values(ii-1)%theta**2)*(ue(ii-1)**6) +  dth2ue6)/(ue(ii)**6))
    end if

    lambda = bl_values(ii)%theta**2.d0*dueds(ii)*Re
    
    ! Check for laminar separation
    if (lambda .lt. -0.0842d0) then
        transorlamsep = 2
        itrans = ii
        exit loop1
    end if
    
    bl_values(ii)%H = fH(lambda)
    L               = fL(lambda)
    
    bl_values(ii)%cf = 2.d0*L/(Re*bl_values(ii)%theta)
    if (ii .gt. 1) then
        bl_values(ii)%cf = bl_values(ii)%cf/ue(ii)
    end if
    
    if ((ii .eq. 1) .and. (dueds(1) .lt. 0.d0)) then
        transorlamsep = 1
        itrans = ii
        exit loop1
    end if
    
    ! Check for transition
    rex = Re*s(ii)*ue(ii)
    ret = Re*bl_values(ii)%theta*ue(ii)
    retmax = 1.174d0*(1.d0 + 22400.d0/rex)*rex**0.46d0
    
    ! Check if rtheta0 is reached
    if ((ii .gt. 2) .and. (scrit .eq. 0.d0)) then
        if (ret .ge. rtheta0(bl_values(ii)%H)) then
            scrit = (s(ii) - s(ii-1))/(rtheta0(bl_values(ii)%H) - rtheta0(bl_values(ii-1)%H))*(ret - rtheta0(bl_values(ii-1)%H)) + s(ii-1)
        end if
    end if
    
    ! Integrate n-factor
    if (scrit .gt. 0.d0) then
        if (bl_values(ii-1)%nfac .eq. 0.d0) then
            bl_values(ii)%nfac = bl_values(ii-1)%nfac + dnds(bl_values(ii)%H,bl_values(ii)%theta)*(s(ii) - scrit)/2.d0
        else
            bl_values(ii)%nfac = bl_values(ii-1)%nfac + (dnds(bl_values(ii)%H,bl_values(ii)%theta) + dnds(bl_values(ii-1)%H,bl_values(ii-1)%theta))/2.d0*(s(ii) - s(ii-1))
        end if
    end if
    
    if (trans_model .eq. 1) then
        if (ret .gt. retmax) then
            transorlamsep = 1
            itrans = ii
            exit loop1
        end if
    else if (trans_model .eq. 2) then
        if (bl_values(ii)%nfac .ge. n_crit) then
            transorlamsep = 1
            itrans = ii
            exit loop1
        end if
    else if (trans_model .eq. 3) then
        if (nodes(ii,1) .ge. x_trans) then
            transorlamsep = 1
            itrans = ii
            exit loop1
        end if
    end if
    
    if (bl_values(ii)%theta /= bl_values(ii)%theta) then
        transorlamsep = 1
        itrans = ii
        exit loop1
    end if
    
    ii = ii + 1
    
    ! Check, if end of airfoil surface is reached
    if (ii .gt. nblp) then
        endofsurf = 1
        itrans    = nblp
        write(*,*) 'itrans end of surface reached'
    end if
end do loop1

! Transition
transloc = 100.d0*nodes(itrans,1)
if (transloc .lt. 0.d0) transloc = 0.1d-3

isep = -1
if (itrans .lt. nblp) then
    ! Turbulent boundary layer
    ! According to HEAD
    ii = itrans

    if (ii .gt. nblp) then
        endofsurf = 1
    end if
    
    H_guess = 1.35d0
!    H_guess = 1.6d0

!    ! Approximation of initial H-value for turbulent boundary
!    ! layer calculation
!    rtheta = Re*bl_values(ii)%theta*ue(ii)
!    if (Re*bl_values(ii)%theta*ue(ii) .lt. 5.d4) then
!        H_delta = 0.821d0 + 0.114d0*log10(rtheta)
!    else
!        H_delta = 1.375d0
!    end if
!    H_guess = bl_values(ii-1)%H - H_delta
    
    y(2) = H1ofH(H_guess)
!    y(2) = H1ofH_houwink(H_guess)
    y(1) = bl_values(ii-1)%theta

    do while ((endofsurf .eq. 0) .and. (turbseparation .eq. 0))

        call runge(s(ii)-s(ii-1),y,Re,ue(ii-1),dueds(ii-1),ue(ii),dueds(ii),temp)
        y(:) = temp(:)
        bl_values(ii)%theta = y(1)
        bl_values(ii)%H = HofH1(y(2))
        
        rtheta = Re*bl_values(ii)%theta*ue(ii)
        bl_values(ii)%cf = cfturb_green(rtheta,bl_values(ii)%H)
!         bl_values(ii)%cf = cfturb(rtheta,bl_values(ii)%H)
        
        if (bl_values(ii)%H /= bl_values(ii)%H) then
            bl_values(ii)%H = 3.d0
        end if

        if (bl_values(ii)%H .ge. 2.4d0) then
            tsep = nodes(ii,1)/nodes(nblp,1)
            turbseparation = 1
            isep = ii
            if (tsep .lt. 0.d0) tsep = 0.d0
            do jj = ii+1, nblp
               bl_values(jj)%H = (bl_values(ii)%H - bl_values(ii-1)%H)/(s(ii) - s(ii-1))*(s(jj) - s(ii-1)) + bl_values(ii-1)%H
               bl_values(jj)%theta = (bl_values(ii)%theta - bl_values(ii-1)%theta)/(s(ii) - s(ii-1))*(s(jj) - s(ii-1)) + bl_values(ii-1)%theta
               bl_values(jj)%cf = (bl_values(ii)%cf - bl_values(ii-1)%cf)/(s(ii) - s(ii-1))*(s(jj) - s(ii-1)) + bl_values(ii-1)%cf
!                bl_values(jj)%H = bl_values(ii)%H
!                bl_values(jj)%theta = bl_values(ii)%theta
!                bl_values(jj)%cf = bl_values(ii)%cf
            end do
            ii = ii - 1
        end if
        
        ii = ii+1

        if (ii .gt. nblp) then
            endofsurf = 1
        end if
    end do
end if

bl_values(:)%H = bl_values(:)%H + 0.2d0*(Ma*ue(:))**2*(bl_values(:)%H + 1.d0)
bl_values(:)%theta = bl_values(:)%theta*(1.d0 + 0.2d0*(Ma*ue(:))**2)**(2.4d0/0.8d0)
bl_values(:)%cf = bl_values(:)%cf*(1.d0 + 0.2d0*(Ma*ue(:))**2)

bl_values(:)%deltas = bl_values(:)%H*bl_values(:)%theta

! Always use point before last point, as the velocity at the last point is not physically based.
if (turbseparation .eq. 0) then
    n_te = nblp-1
else
    n_te = isep-1
    !n_te = nblp - 1
end if

if (turbseparation .eq. 1) then
    if (tsep .eq. 1.d0) tsep = 1.d0 - 1.d-5
    dxsep = nodes(nblp,1) - tsep

    if (isep .ne. nblp) then
        phi_te = atan(-(nodes(isep,2)-nodes(nblp,2))/(nodes(isep,1)-nodes(nblp,1)))
    else
        phi_te = 0.d0
    end if
else
    dxsep  = 0.d0
    phi_te = 0.d0
end if

if (verbose .eqv. .true.) then
    write(*,*)
    write(*,*) 'trans_model:     ', trans_model
    write(*,*) 'theta_te:        ', bl_values(n_te)%theta
    write(*,*) 'H_te:            ', bl_values(n_te)%H
    write(*,*) 'ue_te:           ', ue(n_te)
    write(*,*) 'nblp:            ', nblp
    write(*,*) 'itrans:          ', itrans
    write(*,*) 'transorlamsep:   ', transorlamsep
    write(*,*) 'turbseparation:  ', turbseparation
    write(*,*) 'transloc:        ', transloc
    write(*,*) 'tsep:            ', tsep
    write(*,*) 'dxsep:           ', dxsep
end if

if (bl_values(n_te)%theta /= bl_values(n_te)%theta) then
    write(*,*) 'Boundary layer error occured (theta_te = NaN)'
    error = .true.
    return
end if


contains

    subroutine running_coordinate(n,nodes,s)

    implicit none
    
    integer, intent(in)                          :: n
    double precision, dimension(n,2), intent(in) :: nodes
    double precision, dimension(n), intent(out)  :: s
    integer                                      :: ii

    s(1) = 0.d0
    do ii = 2,n
        s(ii) = s(ii-1) + sqrt((nodes(ii,1) - nodes(ii-1,1))**2 + (nodes(ii,2) - nodes(ii-1,2))**2)
    end do

    end subroutine running_coordinate


    subroutine velocity_gradient(n,ue,s,dueds)

    implicit none

    integer, intent(in)                         :: n
    double precision, dimension(n), intent(in)  :: ue,s
    double precision, dimension(n), intent(out) :: dueds
    integer                                     :: ii
    double precision                            :: v1,v2,v3,x1,x2,x3,fac

    v1 = ue(3)
    v2 = ue(1)
    x1 = s(3)
    x2 = s(1)
    do ii = 1, n
        v3 = v1
        x3 = x1
        v1 = v2
        x1 = x2
        if (ii .lt. nblp) then
            v2 = ue(ii+1)
            x2 = s(ii+1)
        else
            v2 = ue(nblp-2)
            x2 = s(nblp-2)
        end if
        fac = (x3 - x1)/(x2 - x1)
        dueds(ii) = ((v2 - v1)*fac - (v3 - v1)/fac)/(x3 - x2)
    end do
    
    ! If first velocity gradient is negative, recalculate it using
    ! simple different quotient
    if (dueds(1) .lt. 0.d0) then
        dueds(1) = (ue(2) - ue(1))/(s(2) - s(1))
    end if

    end subroutine velocity_gradient


end subroutine solve_bl_moran
