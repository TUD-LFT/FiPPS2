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
!> Benutzerdefinierte Ergebnisausgabe fuer lineare Stabilitaetsanalyse
!
!> @details
!
!> @author Florian Dexl, TU Dresden, WiMi, 28.07.2021
!
!> $Id: sol1_output_hymowi2.f90 434 2021-01-04 16:50:30Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 434 $
!> $Date: 2021-01-04 17:50:30 +0100 (Mo, 04. Jan 2021) $
!
! =================================================================================================
subroutine sol2_output_user(fesim, scloop, eigenvalues, eigenvectors, num_out_eigval, aeroConverged)
!
! use
!
  use fesimulation_typen
  use konstanten
  
  implicit none
  
  type(fe_simulation), intent(in)                         :: fesim
  double precision, dimension(num_out_eigval), intent(in) :: eigenvalues
  double precision, dimension(fesim%internals%dim_dof,num_out_eigval), intent(in) :: eigenvectors
  integer, intent(in)                                     :: num_out_eigval
  integer,intent(in)                                      :: scloop
  logical, intent(in)                                     :: aeroConverged
  
  ! Routine

end subroutine sol2_output_user
