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
! =================================================================================================
!
!> @brief
!> writes computation times to file
!
!> @details
!> writes computation times to text file compTime_Fipps.txt;
!> input values correspond to compiler-function SYSTEM_CLOCK;
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: writeCompTime.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine writeCompTime (compTime,solvTime,solvEigTime,initTime,stiffTime,lvecTime)

! =================================================================================================
! Use
!
use globale_variablen, ONLY : compTimeOutFile
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Include
!
!
! =================================================================================================
!
! Data types
!
! Input
!
double precision, intent(in)                       :: compTime
double precision, intent(in)                       :: solvTime
double precision, intent(in)                       :: stiffTime
double precision, intent(in)                       :: solvEigTime
double precision, intent(in)                       :: initTime
double precision, intent(in)                       :: lvecTime
!
! =================================================================================================
!
! Calculation
!

 OPEN(UNIT=compTimeOutFile,FILE='compTime_Fipps.txt',STATUS='UNKNOWN')
 
 WRITE(compTimeOutFile, '(A164)') ' FiPPS - Computation time [s]     Solver time [s]      Eigsolver time [s]       Init time [s]      Create Stiffness matrix time [s]      Create loadvector time [s]'
 WRITE(compTimeOutFile, '(X,E17.10,16X,E17.10,4X,E17.10,8X,E17.10,2X,E17.10,21X,E17.10)') compTime, solvTime, solvEigTime, initTime, stiffTime, lvecTime

 CLOSE(compTimeOutFile)
end subroutine writeCompTime
