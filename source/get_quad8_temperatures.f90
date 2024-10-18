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
!> Subroutine computes or assignes the temperatures regarding to the specified
!> load case to each element.
!
!> @details
!> Subroutine computes or assignes the temperatures regarding to the specified
!> load case to each element.
!
!> @author Florian Dexl, 22.05.2019
!
!> $Id: get_node_temperatures.f90 394 2019-01-22 07:22:11Z DOM+ahauffe $
!> $Author: DOM+ahauffe $
!> $Revision: 394 $
!> $Date: 2019-01-22 08:22:11 +0100 (Di, 22. Jan 2019) $
!
! =================================================================================================
subroutine get_quad8_temperatures(fesim,lcID, quad8_temperatures)
! =================================================================================================
! use
!
use fesimulation_typen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation)                    :: fesim
integer, intent(in)                    :: lcID !< Load case-ID
!
! Output
double precision, dimension(size(fesim%elemente%quad8s,1)), intent(out)    :: quad8_temperatures    !< Temperatures on each node for this load case
!
! Internal
!
integer                                :: ii, jj
double precision                       :: T_ref

integer                                :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
 quad8_temperatures    = 0.D0
 T_ref = 0.d0
!
! =================================================================================================
!
! Calculation
!
 ! Find reference temperature for current load case
 ! If no reference temperature is given, it is set to 0
 do ii = 1,size(fesim%lasten%quad8temps,1)
   do jj = 1,size(fesim%lasten%loads,1)
     if ((fesim%lasten%loads(jj)%lcid .eq. lcID) .and. (fesim%lasten%loads(jj)%lidi .eq. fesim%lasten%quad8temps(ii)%lid) .and. (fesim%lasten%quad8temps(ii)%eid .eq. 0)) then
       T_ref = fesim%lasten%quad8temps(ii)%temp
     end if
   end do
 end do
 
 ! Loop over each temperature load to get all element temperatures assigned to this load case;
 ! Set temperature at element to (given temperature - reference temperature);
 ! Temperature = 0 is set to all elements where no explicit temperature was assigned
 do ii=1,size(fesim%lasten%quad8temps,1)
   do jj = 1,size(fesim%lasten%loads,1)
     if((fesim%lasten%loads(jj)%lcid .eq. lcID) .and. (fesim%lasten%loads(jj)%lidi .eq. fesim%lasten%quad8temps(ii)%lid) .and. (fesim%lasten%quad8temps(ii)%eid .ne. 0)) then
       ! Assign temperature to element
       quad8_temperatures(fesim%lasten%quad8temps(ii)%eid) = fesim%lasten%quad8temps(ii)%temp - T_ref
     end if
   end do
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'get_quad8_temperatures'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine get_quad8_temperatures
