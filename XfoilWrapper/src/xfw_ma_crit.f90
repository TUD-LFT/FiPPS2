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
SUBROUTINE xfw_ma_crit(n_cp, cp, mach, mach_crit)

implicit none

! Variablendeklaration
intent(in)                                        :: n_cp, cp, mach
intent(out)                                       :: mach_crit
double precision                                  :: mach, mach_crit
double precision                                  :: Ma_eps = 1.d-6
double precision                                  :: cp_min, dMa, func_val!, dfunc_val
double precision, dimension(n_cp)                 :: cp, cp_uncorr
integer                                           :: ii, n_cp, n_iter, iter_max=1000

! Berechne Druckbeiwerte für Ma = 0
! Durch "inverse" Karman-Tsien Korrektur
if (mach .gt. 0.d0) then
    do ii = 1, n_cp
        cp_uncorr(ii) = KarmanTsienInvers(cp(ii),mach)
    end do
else
    cp_uncorr(:) = cp(:)
end if

! Minimaler Druckbeiwert vor Kompressibilitaetskorrektur
cp_min = minval(cp_uncorr)
if (cp_min .ge. 0.d0) then
    write(*,*) 'Minimum pressure coefficient must be less than zero!'
    write(*,*) 'Minimum pressure found: ', cp_min
    STOP
end if

! Berechnung der kritischen Mach-Zahl
mach_crit = 0.5d0
n_iter  = 0
do while (.true.)
    n_iter = n_iter + 1
    func_val = HelpFunc(cp_min, mach_crit)
!     dfunc_val = (HelpFunc(cp_min, mach_crit + Ma_eps) - func_val)/Ma_eps
!     dMa = func_val/dfunc_val
    dMa = func_val/dHelpFunc(cp_min, mach_crit)
    mach_crit = mach_crit - dMa
    if (abs(dMa) .lt. Ma_eps) exit
    if (n_iter .gt. iter_max) then
        write(*,*) 'Calculation of mach_crit failed!'
        STOP
    end if
end do

contains

  double precision FUNCTION KarmanTsien(cp,Ma)
    !
    implicit none

    double precision                               :: cp,Ma

    if (Ma .gt. 0.d0) then
       KarmanTsien = cp/(sqrt(1 - Ma**2.d0) + (Ma**2.d0/(1.d0 + sqrt(1.d0 - Ma**2.d0)))*cp/2.d0)
    else
       KarmanTsien = cp
    end if

  END FUNCTION KarmanTsien

  double precision FUNCTION dKarmanTsien_dMa(cp,Ma)
    !
    implicit none

    double precision                               :: cp,Ma

    if (Ma .gt. 0.d0) then
      dKarmanTsien_dMa = 2.d0*(cp-2.d0)*cp*Ma*(Ma**2-2.d0*(sqrt(1.d0-Ma**2)+1.d0))
      dKarmanTsien_dMa = dKarmanTsien_dMa/(sqrt(1.d0-Ma**2)*((cp-2.d0)*Ma**2+2.d0*(sqrt(1.d0-Ma**2)+1.d0))**2)
    else
      dKarmanTsien_dMa = 1.d0
    end if

  END FUNCTION dKarmanTsien_dMa

  double precision FUNCTION KarmanTsienInvers(cp_KT,Ma)
    !
    implicit none

    double precision                               :: cp_KT,Ma

    if (Ma .gt. 0.d0) then
       KarmanTsienInvers = 2.d0*(Ma**2*cp_KT - sqrt(1.d0-Ma**2)*cp_KT - cp_KT)/(Ma**2*cp_KT - 2.d0*sqrt(1.d0-Ma**2) - 2.d0)
    else
       KarmanTsienInvers = cp_KT
    end if

  END FUNCTION KarmanTsienInvers

  double precision FUNCTION HelpFunc(cpmin,Ma)
    !
    implicit none

    double precision                               :: cpmin,Ma

    HelpFunc = CPcrit(Ma) - KarmanTsien(cpmin,Ma)

  END FUNCTION HelpFunc

  double precision FUNCTION dHelpFunc(cpmin,Ma)
    !
    implicit none

    double precision                               :: cpmin,Ma

    dHelpFunc = dCPcrit_dMa(Ma) - dKarmanTsien_dMa(cpmin,Ma)

  END FUNCTION dHelpFunc

  double precision FUNCTION CPcrit(Ma)
    !
    implicit none

     double precision                               :: Ma
     double precision                               :: kappa=1.4d0

     CPcrit = 1.d0/(0.5d0*kappa*Ma**2)
     CPcrit = CPcrit*(((1.d0+0.5d0*(kappa-1.d0)*Ma**2)/(0.5d0*(kappa+1.d0)))**(kappa/(kappa-1.d0))-1.d0)

  END FUNCTION CPcrit

  double precision FUNCTION dCPcrit_dMa(Ma)
    !
    implicit none

     double precision                               :: Ma
     double precision                               :: kappa=1.4d0

     dCPcrit_dMa = (4.d0*Ma**2-8.d0)*(((kappa-1.d0)*Ma**2+2.d0)/(kappa+1.d0))**(kappa/(kappa-1.d0))+(4.d0*kappa-4.d0)*Ma**2+8.d0
     dCPcrit_dMa = dCPcrit_dMa/(kappa*Ma**3*((kappa-1.d0)*Ma**2+2.d0))

  END FUNCTION dCPcrit_dMa

END SUBROUTINE xfw_ma_crit
