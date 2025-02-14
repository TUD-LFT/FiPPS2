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
!> Subroutine to deallocate all global dynamic arrays
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 06.05.2010
!
!> $Id: free_globals.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine free_globals()

  use globale_variablen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Internal
!
  integer  :: err_code=0

  if (allocated(mpi_quad8_proc_dist) .eq. .true.) deallocate(mpi_quad8_proc_dist)
  if (allocated(mpi_lsolid20_proc_dist) .eq. .true.) deallocate(mpi_lsolid20_proc_dist)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                       'An error occured in subroutine'
   write(*,*)                       'Subroutinenname'
   write(*,'(A,I2)',advance='YES')  ' Errorcode: ', err_code
   write(*,*)                       'exit program '
   stop
   
end if  
end subroutine free_globals
