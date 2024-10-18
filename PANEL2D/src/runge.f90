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
SUBROUTINE runge(dx,y,Re,uei,duedsi,ueip1,duedsip1,ynp1)

use konst_var
use bl_functions

implicit none

! Declaration of variables
intent(in)                                        :: dx,y,Re,uei,duedsi,ueip1,duedsip1
intent(out)                                       :: ynp1
double precision                                  :: dx,Re,uei,duedsi,ueip1,duedsip1
double precision, dimension(2)                    :: y,ynp1

double precision, dimension(2)                    :: yt,yp,ys
double precision                                  :: H1,H,rtheta

double precision, parameter                       :: eps = 1.d-10

integer                                           :: tsep

! Calculation
tsep = 0

! 1st step
yt(:) = y(:)

H1 = yt(2)
H = HofH1(H1)

! if ((H .eq. 3.d0) .or. (H1 .le. 3.d0)) then
if (H1 .le. 3.d0)  then
    tsep = 1
end if

if (tsep .eq. 0) then
    rtheta = Re*uei*yt(1)
    
    yp(1) = -(H + 2.d0)*yt(1)*duedsi/uei + 0.5d0*cfturb_green(rtheta,H)
!     yp(1) = -(H + 2.d0)*yt(1)*duedsi/uei + 0.5d0*cfturb(rtheta,H)
    yp(2) = -H1*(duedsi/uei + yp(1)/yt(1)) + 0.0306d0*(H1 - 3.d0)**(-0.6169d0)/yt(1)
    
    yt(:) = y(:) + dx*yp(:)
    
    ys(:) = y(:) + 0.5d0*dx*yp(:)
    
    ! 2nd step
    H1 = yt(2)
    H  = HofH1(H1)

!     if ((H .eq. 3.d0) .or. (H1 .le. 3.d0)) then
    if ((abs(H1 - 3.d0) .lt. eps) .or. (H1 .lt. 3.d0))  then
        tsep = 1
    end if

    if (tsep .eq. 0) then
        rtheta = Re*ueip1*yt(1)
        
        yp(1) = -(H + 2.d0)*yt(1)*duedsip1/ueip1 + 0.5d0*cfturb_green(rtheta,H)
!         yp(1) = -(H + 2.d0)*yt(1)*duedsip1/ueip1 + 0.5d0*cfturb(rtheta,H)
        yp(2) = -H1*(duedsip1/ueip1 + yp(1)/yt(1)) + 0.0306*(H1 - 3.d0)**(-0.6169d0)/yt(1)

        ynp1(:) = ys(:) + 0.5d0*dx*yp(:)
     end if
end if

if (tsep .eq. 1) then
    ynp1(1) = -2.d0     ! Letting H = 3.d0 in the main routine
    ynp1(2) = 0.d0
end  if
 
END SUBROUTINE runge
