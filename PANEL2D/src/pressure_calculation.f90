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
SUBROUTINE pressure_calculation(ue, cp, p_mano, C_lift, C_mom, Ma_crit)

! Load modules
use konst_var
use netz_variablen

implicit none

! Declaration of variables
double precision, dimension(n_panels), intent(in)  :: ue
double precision, dimension(n_panels), intent(out) :: cp,p_mano
double precision, intent(out)                      :: C_lift, C_mom
double precision, intent(out), optional            :: Ma_crit
double precision                                   :: Ma_eps = 1.d-6
double precision                                   :: cp_min, dMa, func_val!, dfunc_val
double precision, parameter                        :: damp=7.d-1
double precision, dimension(2), parameter          :: mom_ref=(/0.25d0,0.d0/)
integer                                            :: ii, n_iter, iter_max=10000

! Calculate pressure coefficients
cp(:) = 1.d0 - ue(:)**2.d0

! Minimal pressure coefficient before compressibility correction
cp_min = minval(cp)
if (cp_min .ge. 0.d0) then
    write(*,*) 'Minimum pressure coefficient must be less than zero!'
    write(*,*) 'Minimum pressure found: ', cp_min
    stop
end if

! Calculate lift coefficient
C_lift = 0.d0
do ii = 1, n_panels
    C_lift = C_lift - cp(ii)*lengths(ii)*dot_product(n_vec(ii,:),(/-sin(aoa),cos(aoa)/))
end do

! Calculate moment coefficient
C_mom = 0.d0
do ii = 1, n_panels
    C_mom = C_mom + cp(ii)*lengths(ii)*dot_product(n_vec(ii,:),(/-coll(ii,2)+mom_ref(2),coll(ii,1)-mom_ref(1)/))
end do

! Prandtl-Glauert correction of pressure coefficients, if Mach > 0
do ii = 1, n_panels
    cp(ii) = PrandtlGlauert(cp(ii),Ma)
end do

! Prandtl-Glauert correction of the lift coefficient, if Mach > 0
C_lift = PrandtlGlauert(C_lift,Ma)

! Prandtl-Glauert correction of the moment coefficient, if Mach > 0
C_mom  = PrandtlGlauert(C_mom,Ma)

! Calculate manometer pressure values
p_mano(:) = cp(:)*q_infty

! Calculation of critical Mach number
if (present(Ma_crit)) then
  Ma_crit = 0.5d0
  n_iter  = 0
  do while (.true.)
    n_iter = n_iter + 1
    func_val = HelpFunc(cp_min, Ma_crit)
    !dfunc_val = (HelpFunc(cp_min, Ma_crit + Ma_eps) - func_val)/Ma_eps
    !dMa = func_val/dfunc_val
    dMa = func_val/dHelpFunc(cp_min, Ma_crit)
    Ma_crit = Ma_crit - dMa*damp
    if (abs(dMa) .lt. Ma_eps) exit
    if (n_iter .gt. iter_max) then
        write(*,*) 'Calculation of Ma_crit failed!'
        STOP
    end if
  end do
  if (Ma_crit .lt. 0.d0) Ma_crit = 0.d0
end if

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
  
  double precision FUNCTION PrandtlGlauert(cp,Ma)
    !
    implicit none

    double precision                               :: cp,Ma

     if (Ma .gt. 0.d0) then
        PrandtlGlauert = cp/sqrt(1.d0 - Ma**2)
     else
        PrandtlGlauert = cp
    end if

  END FUNCTION PrandtlGlauert
  
  double precision FUNCTION dPrandtlGlauert_dMa(cp,Ma)
    !
    implicit none

    double precision                               :: cp,Ma

     if (Ma .gt. 0.d0) then
        dPrandtlGlauert_dMa = cp*Ma/(sqrt(1.d0 - Ma**2))**1.5d0
     else
        dPrandtlGlauert_dMa = 1.d0
    end if

  END FUNCTION dPrandtlGlauert_dMa

  double precision FUNCTION HelpFunc(cpmin,Ma)
    !
    implicit none

    double precision                               :: cpmin,Ma

    HelpFunc = CPcrit(Ma) - PrandtlGlauert(cpmin,Ma)

  END FUNCTION HelpFunc

  double precision FUNCTION dHelpFunc(cpmin,Ma)
    !
    implicit none

    double precision                               :: cpmin,Ma

    dHelpFunc = dCPcrit_dMa(Ma) - dPrandtlGlauert_dMa(cpmin,Ma)

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

END SUBROUTINE pressure_calculation
