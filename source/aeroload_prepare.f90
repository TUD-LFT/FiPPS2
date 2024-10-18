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
!> Preperation of aeroloads
!
!> @details
!> Preperation of aeroloads. Determine the internal aeroload IDs and corresponding load factors
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 11.09.2018
!
!> $Id: aeroload_prepare.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine aeroload_prepare(fesim)
! =================================================================================================
!
! Use
!
  use fesimulation_typen
  
  type(fe_simulation) :: fesim
  
  logical :: found
  
  integer :: hh,jj
  integer :: err_code 
  
  err_code = 0
  
!--------------------------------------------------------------------------------------------------
!

fesim%lasten%subcases(:)%aeroloadID  = 0
fesim%lasten%subcases(:)%aeroloadFac = 1.d0

if ((fesim%is_aeroload2d .eq. .true.) .and. (fesim%is_aeroload3d .eq. .true.)) then
  write(*,*) "aeroload2d aeroload3d may not"
  write(*,*) "be used in combination."
  err_code=1
  goto 9999
end if

! Aeroload2D
if (fesim%is_aeroload2d .eq. .true.) then
  do scloop = 1,fesim%num_subcases
    found = .false.
    do hh = 1,size(fesim%lasten%loads,1)
      if (fesim%lasten%loads(hh)%lcid .eq. fesim%lasten%subcases(scloop)%loadid) then
        do jj=1,size(fesim%lasten%aeroload2ds,1)
          if (fesim%lasten%loads(hh)%lidi == fesim%lasten%aeroload2ds(jj)%lid) then
          
            if (found .eqv. .true.) then
              write(*,*) 'Je Subcase darf nur eine aerodynamische Last aeroload definiert sein!'
              err_code=1
              goto 9999
            end if
            
            fesim%lasten%subcases(scloop)%aeroloadID  = jj
            fesim%lasten%subcases(scloop)%aeroloadFac = fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
            
            found = .true.
            
          end if
        end do
      end if
    end do
  end do
end if

! Aeroload3D
if (fesim%is_aeroload3d .eq. .true.) then
  do scloop = 1,fesim%num_subcases
    found = .false.
    do hh = 1,size(fesim%lasten%loads,1)
      if (fesim%lasten%loads(hh)%lcid .eq. fesim%lasten%subcases(scloop)%loadid) then
        do jj=1,size(fesim%lasten%aeroload3ds,1)
          if (fesim%lasten%loads(hh)%lidi == fesim%lasten%aeroload3ds(jj)%lid) then
          
            if (found .eqv. .true.) then
              write(*,*) 'Je Subcase darf nur eine aerodynamische Last aeroload definiert sein!'
              err_code=1
              goto 9999
            end if
            
            fesim%lasten%subcases(scloop)%aeroloadID  = jj
            fesim%lasten%subcases(scloop)%aeroloadFac = fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
            
            found = .true.
            
          end if
        end do
      end if
    end do
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
    write(*,*)                      'aeroload_prepare'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine aeroload_prepare
