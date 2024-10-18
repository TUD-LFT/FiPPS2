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
SUBROUTINE isa(H, T, p, rho, a, nu)

implicit none

! Variablendeklaration
intent(in)                                        :: H
intent(out)                                       :: T, p, rho, a, nu

double precision                                  :: H, T, p, rho, a, nu

double precision                                  :: H_A, gamma, T_A, p_A, rho_A
double precision, parameter                       :: g = 9.80665d0
double precision, parameter                       :: R = 2.8705d2
double precision, parameter                       :: k = 1.4d0
double precision, parameter                       :: nu_0 = 1.7894d-5
double precision, parameter                       :: T_0 = 2.7311d2
double precision, parameter                       :: S = 1.1056d2

! Umgebungsbedingungen nach internationaler Standardatmosphaere
if ((H .gt. -5.d3) .and. (H .le. 1.1d4)) then
    H_A   = 0.d0
    gamma = -6.5d-3
    T_A   = 2.8815d2
    p_A   = 1.01325d5
    rho_A = 1.225d0
else if ((H .gt. 1.1d4) .and. (H .le. 2.0d4)) then
    H_A   = 1.1d4
    gamma = 0.d0
    T_A   = 2.1665d2
    p_A   = 2.262d4
    rho_A = 3.637d-1
else if ((H .gt. 2.0d4) .and. (H .le. 3.2d4)) then
    H_A   = 2.0d4
    gamma = 1.d-3
    T_A   = 2.1665d2
    p_A   = 5.469d3
    rho_A = 8.79d-2
else
    write(*,*)
    write(*,*) '  FEHLER: Standardatmosphaere fuer H = ', H, ' m nicht definiert!'
    write(*,*)
    STOP
end if

if (((H .GT. -5.d3) .AND. (H .LE. 1.1d4)) .OR. ((H .GT. 2.d4) .AND. (H .LE. 3.2d4))) then
  T = T_A + gamma*(H - H_A)
  p = p_A*(1.d0 + gamma/T_A*(H - H_A))**(-1.d0*g/(R*gamma))
  rho = rho_A*(1.d0 + gamma/T_A*(H - H_A))**(-1.d0*g/(R*gamma) - 1.d0)
else
  T = T_A
  p = p_A*exp(-1.d0*g*(H-H_A)/(R*T_A))
  rho = rho_A*exp(-1.d0*g*(H-H_A)/(R*T_A))
end if

! Schallgeschwindigkeit berechnen
a = sqrt(k*R*T)

! Berechnung der dynamsichen Viskositaet nach Sunderland-Gesetz
! mit drei Koeffizienten
nu = nu_0*(T/T_0)**1.5d0*(T_0 + S)/(T + S)

END SUBROUTINE isa