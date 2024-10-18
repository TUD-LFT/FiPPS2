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
!> Benutzerdefinierte Ergebnisausgabe fuer statische Loesung
!
!> @details
!
!> @author Florian Dexl, TU Dresden, WiMi, 27.07.2021
!
!> $Id: sol1_output_hymowi2.f90 434 2021-01-04 16:50:30Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 434 $
!> $Date: 2021-01-04 17:50:30 +0100 (Mo, 04. Jan 2021) $
!
! =================================================================================================
subroutine sol1_output_user(fesim,Utot,Fout,scloop,cdi,aeroConverged)
!
! use
!
  use fesimulation_typen
  use vtk_variablen
  use failure_criteria
  
  implicit none
  
  type(fe_simulation), intent(in)                         :: fesim
  double precision, dimension(fesim%num_dof), intent(in)  :: Utot, Fout
  integer,intent(in)                                      :: scloop
  double precision, intent(in)                            :: cdi
  logical, intent(in)                                     :: aeroConverged
  
  integer                                                 :: err_code=0
  
  ! Routine
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_user'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if  
  
end subroutine sol1_output_user
