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
SUBROUTINE boundary_layer(verbose,n_points,points,ue,re,C_drag,DeltaCL,FactorCM,dxsep,transloc_up,transloc_lo,bl_values,error,C_lift)

! Load modules
use konst_var
use functions

implicit none

! Declaration of variables
intent(inout)                                     :: ue
intent(in)                                        :: verbose,n_points,points,re,C_lift
intent(out)                                       :: C_drag,DeltaCL,FactorCM,dxsep,transloc_up,transloc_lo,bl_values,error
integer                                           :: n_points
double precision                                  :: re,C_drag,DeltaCL,FactorCM,dxsep,C_lift
double precision, dimension(n_points)             :: ue
type(bl_point), dimension(n_points)               :: bl_values
double precision, dimension(n_points,2)           :: points
logical                                           :: verbose

double precision, dimension(:,:), allocatable     :: zup,zlo
double precision, dimension(:), allocatable       :: ueup,uelo
type(bl_point), dimension(:), allocatable         :: bl_up, bl_lo
double precision, dimension(2)                    :: zstag
double precision                                  :: tsep_up, tsep_lo
double precision                                  :: dxsep_up, dxsep_lo
double precision                                  :: transloc_up, transloc_lo
double precision                                  :: DeltaCD
double precision                                  :: ue_stag
double precision                                  :: phi_te_lo, phi_te_up, DeltaCL_up, DeltaCL_lo
integer                                           :: n_te_up, n_te_lo

integer                                           :: ispu,ispl,nup,nlo

logical                                           :: error

! Find stagnation point next to airfoil's leading edge
call stagnation_point(ue, ispu, ispl)

! Division of boundary layer in upper and lower airfoil side
nlo = ispl
nup = n_points+1-ispu
if (ispu .ne. ispl) then
    nlo = nlo + 1
    nup = nup + 1
    ue_stag = 0.d0
    zstag(:) = (points(ispu,:) - points(ispl,:))*(ue_stag - ue(ispl))/(ue(ispu) - ue(ispl)) + points(ispl,:)
end if

allocate(zup(nup,2), zlo(nlo,2))
allocate(ueup(nup),  uelo(nlo))
allocate(bl_up(nup), bl_lo(nlo))

if (ispu .ne. ispl) then
    zlo(2:,:) = points(ispl:1:-1,:)
    zup(2:,:) = points(ispu:n_points,:)
    uelo(2:) = ue(ispl:1:-1)
    ueup(2:) = ue(ispu:n_points)
    zlo(1,:) = zstag(:)
    zup(1,:) = zstag(:)
    uelo(1) = ue_stag
    ueup(1) = ue_stag
else
    zlo(:,:) = points(ispu:1:-1,:)
    zup(:,:) = points(ispu:n_points,:)
    uelo(:) = ue(ispl:1:-1)
    ueup(:) = ue(ispl:n_points)
end if
uelo(:) = -1.d0*uelo(:)

! Boundary layer calculation at airfoil's upper side
call solve_bl_moran(verbose,re,zup,nup,ueup,dxsep_up,transloc_up,n_te_up,bl_up,phi_te_up,error)
if (error .eqv. .true.) then
    return
end if

! Boundary layer calculation at airfoil's lower side
call solve_bl_moran(verbose,re,zlo,nlo,uelo,dxsep_lo,transloc_lo,n_te_lo,bl_lo,phi_te_lo,error)
if (error .eqv. .true.) then
    return
end if

!! Correction of drag and lift coefficient according to JavaFOIL
!if (aoa .ge. 0.d0) then
!    DeltaCL = -0.2d0*dxsep_up!*2.d1
!    DeltaCD = (sin(aoa)**2.d0 + abs(0.025d0*cos(aoa)))*(dxsep_up)**2.d0
!else
!    DeltaCL = -0.2d0*dxsep_lo!*2.d1
!    DeltaCD = (sin(aoa)**2.d0 + abs(0.025d0*cos(aoa)))*(dxsep_lo)**2.d0
!end if

! Modified correction of drag and lift coefficient based on JavaFOIL
! without distinguishing between positive and negative angle of attack
!DeltaCL = -0.2d0*(dxsep_up + dxsep_lo)!*2.d0!*1.d2
!DeltaCL = -sin(aoa)*(dxsep_up + dxsep_lo)

if (C_lift .gt. 0.d0) then
  DeltaCL_up = -(aoa + phi_te_up)*PI*dxsep_up
  if (DeltaCL_up .gt. 0.d0) then
    DeltaCL_up = -abs(sin(aoa))*dxsep_up
  end if
  DeltaCL_lo = -abs(sin(aoa))*dxsep_lo
else
  DeltaCL_lo = (aoa + phi_te_lo)*PI*dxsep_lo
  if (DeltaCL_lo .gt. 0.d0) then
    DeltaCL_lo = -abs(sin(aoa))*dxsep_lo
  end if
  DeltaCL_up = -abs(sin(aoa))*dxsep_up
end if
DeltaCL = DeltaCL_lo + DeltaCL_up
FactorCM = 0.9d0*(1.d0-dxsep_up)**2.d0*(1.d0-dxsep_lo)**2.d0
DeltaCD = (sin(aoa)**2.d0 + abs(0.025d0*cos(aoa)))*((dxsep_up)**2.d0 + (dxsep_lo)**2.d0)

! Calculation of drag coefficient according to SQUIRE-YOUNG
call drag_squire_young(bl_up(n_te_up)%theta,bl_up(n_te_up)%H,ueup(n_te_up),bl_lo(n_te_lo)%theta,bl_lo(n_te_lo)%H,uelo(n_te_lo),C_drag)

if (C_drag .le. 0.d0) then
    write(*,*) 'Negative or zero drag!'
    STOP
end if

C_drag = abs(C_drag + DeltaCD)

if (verbose .eqv. .true.) then
    if (dxsep_up .gt. 0.d0) write(*,*) 'Upper side: Turbulent separation occured at ', (1.d0 - dxsep_up)
    if (dxsep_lo .gt. 0.d0) write(*,*) 'Lower side: Turbulent separation occured at ', (1.d0 - dxsep_lo)
    write(*,*)
    write(*,*) 'DeltaCD:           ', DeltaCD
    write(*,*) 'DeltaCL:           ', DeltaCL
    write(*,*) 'FactorCM:          ', FactorCM
    write(*,*) 'ue_te_up:          ', ueup(n_te_up)
    write(*,*) 'th_up:             ', bl_up(n_te_up)%theta
    write(*,*) 'H_up:              ', bl_up(n_te_up)%H
    write(*,*) 'ue_te_lo:          ', uelo(n_te_lo)
    write(*,*) 'th_lo:             ', bl_lo(n_te_lo)%theta
    write(*,*) 'H_lo:              ', bl_lo(n_te_lo)%H
    write(*,*)
    write(*,*) 'C_drag:            ', C_drag
end if

if (ispl .ne. ispu) then
    bl_values(1:(nlo-1))      = bl_lo(nlo:2:-1)
    bl_values(nlo:n_points)   = bl_up(2:nup)
    ue(ispl:1:-1)             = -1.d0*uelo(2:)
    ue(ispu:n_points)         = ueup(2:)
else
    bl_values(1:nlo)          = bl_lo(nlo:1:-1)
    bl_values(nlo+1:n_points) = bl_up(2:nup)
    ue(ispl:1:-1)             = -1.d0*uelo(:)
    ue(ispu:n_points)         = ueup(:)
end if

dxsep = maxval((/dxsep_up,dxsep_lo/))

deallocate(zup,zlo,ueup,uelo,bl_up,bl_lo)

END SUBROUTINE boundary_layer
