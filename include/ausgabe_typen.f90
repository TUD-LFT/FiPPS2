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
!
!> $Id: ausgabe_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
module ausgabe_typen

  implicit none
  
type ausgabe_type
    logical                         :: outputVTK
    logical                         :: outputUser
    logical                         :: outputShort
    logical                         :: outputKoopt
    logical                         :: outputAdviLa
    logical                         :: outputOptitube
    logical                         :: outputGlawi
    logical                         :: outputHyMoWi
    logical                         :: outputElemCoord
    logical                         :: outputBoundCond
    logical                         :: outputApamePressures
  
    double precision                :: xmax =  HUGE(0.d0)
    double precision                :: xmin = -HUGE(0.d0)
    double precision                :: ymax =  HUGE(0.d0)
    double precision                :: ymin = -HUGE(0.d0)
    double precision                :: zmax =  HUGE(0.d0)
    double precision                :: zmin = -HUGE(0.d0)
end type ausgabe_type

end module ausgabe_typen
