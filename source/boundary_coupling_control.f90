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
!> control subroutine for inserting couplings
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 10.10.2012
!
!> $Id: boundary_coupling_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine boundary_coupling_control (fesim, fixedU)

use konstanten
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
type(fe_simulation)                     :: fesim
!
! Output
!
logical, dimension(fesim%num_dof),intent(out) :: fixedU
!
! inner
!
integer                                 :: ii
integer                                 :: err_code=0
integer                                 :: oldCpSet

! Loop over Couplings
if (fesim%is_coupling .EQ. .true.) then
  oldCPSet = fesim%randbedingungen%couplings(1)%cpsid
  do ii = 2, size(fesim%randbedingungen%couplings,1)
    if (fesim%randbedingungen%couplings(ii)%cpsid .eq. oldCPSet) then
      if (fixedU((fesim%randbedingungen%couplings(ii)%nid-1)*ndof+fesim%randbedingungen%couplings(ii)%dof) == .TRUE.) then
         
         write(*,*) 'You try to set a coupling on a dof'
         write(*,*) 'that already is locked via a spc constraints'
         err_code=1
         goto 9999

      else
         
         fixedU((fesim%randbedingungen%couplings(ii)%nid-1)*ndof+fesim%randbedingungen%couplings(ii)%dof) = .TRUE.
         
      end if
    else
      oldCPSet = fesim%randbedingungen%couplings(ii)%cpsid
    end if
  end do
end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'boundary_coupling_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine boundary_coupling_control
