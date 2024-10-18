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
!> calculates the missing value of mat1-card
!
!> @details
!> calculates the missing value of mat1-card if only 2 of the 3 values E,G and
!> nu are specified in the mat1-card
!> if all 3 values are defined on the mat1-cardf these values are returned
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 24.06.2010
!
!> $Id: mat1_calc_missing_value.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine mat1_calc_missing_value(fesim,ii,Emat,Gmat,nu)

use fesimulation_typen
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
type(fe_simulation), intent(in) :: fesim
integer, intent(in)             :: ii        !< mat1-card index values shall be calculated for
!
! Output
!
double precision, intent(out)   :: Emat      !< youngs modulus
double precision, intent(out)   :: Gmat      !< shear modulus
double precision, intent(out)   :: nu        !< poisson ration
!
! Input + Output
!

!
! Internal
!
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
  if (fesim%materialien%mat1s(ii)%ym == 0.D0 .or. fesim%materialien%mat1s(ii)%sm == 0.D0 .or. fesim%materialien%mat1s(ii)%nu == 0.D0) then
    
    if (fesim%materialien%mat1s(ii)%ym /= 0.D0 .and. fesim%materialien%mat1s(ii)%sm /= 0.D0) then
       Emat = fesim%materialien%mat1s(ii)%ym
       Gmat = fesim%materialien%mat1s(ii)%sm
       nu   = Emat/(2.D0*Gmat)-1.D0
    else if (fesim%materialien%mat1s(ii)%ym /= 0.D0 .and. fesim%materialien%mat1s(ii)%nu /= 0.D0) then
       Emat = fesim%materialien%mat1s(ii)%ym
       nu   = fesim%materialien%mat1s(ii)%nu
       Gmat = Emat/(2.D0*(1.D0+nu))
    else if (fesim%materialien%mat1s(ii)%sm /= 0.D0 .and. fesim%materialien%mat1s(ii)%nu /= 0.D0) then
       Gmat = fesim%materialien%mat1s(ii)%sm
       nu   = fesim%materialien%mat1s(ii)%nu
       Emat = 2.D0*(1.D0+nu)*Gmat
    else
       write(*,*) 'Too few input data on mat1-Cards'
       err_code=1
       goto 9999
    end if
    
  else
  
    Emat = fesim%materialien%mat1s(ii)%ym
    Gmat = fesim%materialien%mat1s(ii)%sm
    nu   = fesim%materialien%mat1s(ii)%nu
    
  end if
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'mat1_calc_missing_value'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine mat1_calc_missing_value
