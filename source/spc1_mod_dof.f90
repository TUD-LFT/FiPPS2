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
!> subroutine for calculating the vector of dofs to be locked in a spc1 card
!
!> @details
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit 08.12.2009
!
!> $Id: spc1_mod_dof.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine spc1_mod_dof (count_dof,dof,mod_dof)

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
integer,intent(in)                              :: count_dof, dof
!
! Output
!
integer,dimension(count_dof),intent(out)        :: mod_dof
!
! Input + Output
!

!
! inner
!
integer                                         :: bccomb, size_dof

integer                                         :: ii
integer                                         :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
!

bccomb = dof

do ii=1,count_dof

  size_dof = 10**(count_dof-ii)
  
  mod_dof(ii) = INT(bccomb/DBLE(size_dof))
  
  bccomb = bccomb - mod_dof(ii) * size_dof
  
end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'spc1_mod_dof'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine spc1_mod_dof
