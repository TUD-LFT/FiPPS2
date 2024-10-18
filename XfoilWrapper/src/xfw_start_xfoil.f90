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
SUBROUTINE xfw_start_xfoil(fname_xb)

  use xfw_konst_var, only: start_xfoil

  implicit none
  
  character(len=256), intent(in)                        :: fname_xb

!--------------------------------------------------------------------------------------------------------------------------
!***************************************************************************************************
!               Starte Xfoil
!**************************************************************************************************

  CALL SYSTEM (trim(start_xfoil) // ' < ' // trim(fname_xb) // ' > /dev/null') !// ' > xf_logf.txt')

END SUBROUTINE xfw_start_xfoil
