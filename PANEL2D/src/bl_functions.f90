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
MODULE bl_functions

use konst_var

implicit none

save

contains

double precision FUNCTION fH(lambda)
    !
    implicit none

    double precision                                :: lambda

    if (lambda .lt. 0.d0) then
        fH = 2.088d0 + 0.0731d0/(lambda + 0.14d0)
    else
        fH = 2.61d0 - lambda*(3.75d0 - 5.24d0*lambda)
    end if
END FUNCTION fH


double precision FUNCTION fL(lambda)
    !
    implicit none

    double precision                                :: lambda

    if (lambda .lt. 0.d0) then
        fL = 0.22d0 + 1.402d0*lambda + 0.018d0*lambda/(lambda + 0.107d0)
    else
        fL = 0.22d0 + lambda*(1.57d0 - 1.8d0*lambda)
    end if
END FUNCTION fL


double precision FUNCTION HofH1(H1)
    !
    implicit none

    double precision                                :: H1

    if (H1 .lt. 3.3d0) then
        HofH1 = 3.d0
    else if (H1 .lt. 5.3d0) then
        HofH1 = 0.6778d0 + 1.1536d0*(H1 - 3.3d0)**(-0.326d0)
    else
        HofH1 = 1.1d0 + 0.86d0*(H1 - 3.3d0)**(-0.777d0)
    end if
END FUNCTION HofH1


double precision FUNCTION H1ofH(H)
    !
    implicit none

    double precision                                :: H

!     if (H .lt. 1.1d0) then
!         write(*,*) 'H < 1.1 ! -> H1 = 16'
!         H1ofH = 16.d0
!     else if (H .le. 1.6d0) then
    if (H .le. 1.6d0) then
        H1ofH = 3.3d0 + 0.8234d0*(H - 1.1d0)**(-1.287d0)
    else
        H1ofH = 3.3d0 + 1.5501d0*(H - 0.6778d0)**(-3.064d0)
    end if
END FUNCTION H1ofH


double precision FUNCTION H1ofH_houwink(H)
    !
    implicit none
    
    double precision                                :: H,ht
    
    ! Closure of H1(H) according to Houwink and Veldman
    ! See: Veldman, A.E.P.: A simple interaction law for viscous-
    !      inviscid interaction, J. Eng. Math. 65; S. 367-383; 2009
    !      (Eq. (4))
    ht = min(H, ((H - 2.732d0)/2.d0 + 2.732d0))
    if (H .le. 4.d0) then
        H1ofH_houwink = (ht*(ht + 2.d0))/(2.d0*(ht - 1.d0))
    else
        H1ofH_houwink = 1.75d0 + 5.52273*ht/(ht + 5.818181d0)
    end if
END FUNCTION H1ofH_houwink


double precision FUNCTION cfturb(rtheta,H)
    !
    implicit none
    
    double precision                                :: rtheta,H

    ! Ludwieg-Tillman Skin Friction Formula
    cfturb = 0.246d0*(10.d0**(-0.678d0*H))*rtheta**(-0.268d0)
END FUNCTION cfturb


double precision FUNCTION cfturb_green(rtheta,H)
    !
    implicit none
    
    double precision                                :: rtheta,H
    double precision                                :: cf0, h0
    
    ! Green's modification of Ludwieg-Tillman relationship
    ! See: Veldman, A.E.P.: A simple interaction law for viscous-
    !      inviscid interaction, J. Eng. Math. 65; S. 367-383; 2009
    !      (Eq. (5))
    ! Compared to Ludwieg-Tillman formula, this correlation can
    ! provide negative values for c_f. According to Moran, c_f <= 0
    ! indicates turbulent separation
    cf0          = 0.01013d0/(log10(rtheta) - 1.02d0) - 0.00075d0
    h0           = 1.d0 - 6.55d0*sqrt(0.5d0*cf0)
    cfturb_green = cf0*(0.9d0/(H*h0 - 0.4d0) - 0.5d0)
END FUNCTION cfturb_green


double precision FUNCTION rtheta0(H)
    !
    implicit none

    double precision                                :: H
    
    rtheta0 = 1.415d0/(H - 1.d0) - 0.489d0
    rtheta0 = rtheta0*tanh(20.d0/(H - 1.d0) - 12.9d0)
    rtheta0 = rtheta0 + 3.295d0/(H - 1.d0) + 0.44d0
    rtheta0 = 10.d0**rtheta0
END FUNCTION rtheta0


double precision FUNCTION dnds(H,theta)
    !
    implicit none
    
    double precision                                :: H,theta

    dnds = dndRe(H)*(mfun(H) + 1.d0)/2.d0*lfun(H)/theta
END FUNCTION dnds
    

double precision FUNCTION dndRe(H)
    !
    implicit none
    
    double precision                                :: H

    dndRe = 0.01d0*sqrt((2.4d0*H - 3.7d0 + 2.5d0*tanh(1.5d0*H - 4.65d0))**2.d0 + 0.25d0)
END FUNCTION dndRe

    
double precision FUNCTION lfun(H)
    !
    implicit none
    
    double precision                                :: H

    lfun = (6.54d0*H - 14.07d0)/H**2.d0
END FUNCTION lfun


double precision FUNCTION mfun(H)
    !
    implicit none
    
    double precision                                :: H

    mfun = (0.058d0*(H - 4.d0)**2.d0/(H-1.d0) - 0.068d0)/lfun(H)
END FUNCTION mfun


double precision FUNCTION f1fun(H)
    !
    implicit none

    double precision                                :: H
    
    if (H .lt. 7.4d0) then
        f1fun = -0.067d0 + 0.01977d0*(7.4d0 - H)**2.d0/(H - 1.d0)
    else
        f1fun = -0.067d0 + 0.022d0*(1.d0 - 1.4d0/(H - 6.d0))**2.d0
    end if
END FUNCTION f1fun


double precision FUNCTION df1fun(H)
    !
    implicit none

    double precision                                :: H
    
    if (H .lt. 7.4d0) then
        df1fun = (0.01977*H**2.d0 - 0.03954*H - 0.790009)/(H - 1.d0)**2.d0
    else
        df1fun = (0.0616*H - 0.45584)/(H - 6.d0)**3.d0
    end if
END FUNCTION df1fun


double precision FUNCTION f2fun(H)
    !
    implicit none

    double precision                                :: H

    if (H .lt. 4.d0) then
        f2fun = 0.207d0 + 0.00205d0*(4.d0 - H)**5.5d0
    else
        f2fun = 0.207d0 - 0.003d0*(4.d0 - H)**2.d0/(1.d0 + 0.02d0*(H - 4.d0)**2.d0)
    end if
END FUNCTION f2fun


double precision FUNCTION df2fun(H)
    !
    implicit none

    double precision                                :: H

    if (H .lt. 4.d0) then
        df2fun = -0.011275*(4.d0 - H)**4.5d0
    else
        df2fun = -(15.d0*H - 60.d0)/(H**2.d0 - 8.d0*H + 66.d0)**2.d0
    end if
END FUNCTION df2fun


double precision FUNCTION f3fun(H)
    !
    implicit none

    double precision                                :: H

    f3fun = f2fun(H) - f1fun(H)
END FUNCTION f3fun


double precision FUNCTION df3fun(H)
    !
    implicit none

    double precision                                :: H

    df3fun = df2fun(H) - df1fun(H)
END FUNCTION df3fun


double precision FUNCTION Hs(H)
    !
    implicit none

    double precision                                :: H

    if (H .lt. 4.d0) then
        Hs = 1.515d0 + 0.076d0*(4.d0 - H)**2.d0/H
    else
        Hs = 1.515d0 + 0.04d0*(4.d0 - H)**2.d0/H
    end if
end FUNCTION Hs


double precision FUNCTION gfun(H)
    !
    implicit none

    double precision                                :: H

    if (H .lt. 4.d0) then
        gfun = (H**2.d0 - 16.d0)/(H*(H**2.d0 + 11.9342d0*H + 16.d0))
    else
        gfun = (H**2.d0 - 16.d0)/(H*(H**2.d0 + 29.875d0*H + 16.d0))
    end if
end FUNCTION gfun


END MODULE bl_functions
