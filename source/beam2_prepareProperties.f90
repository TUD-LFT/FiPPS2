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
!
!> @details
!
!> @author 
!
!> $Id: beam2_prepareProperties.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_prepareProperties(fesim)

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
!
! Output
!
!
! Input + Output
!
  type (fe_simulation), intent(inout) :: fesim
!
! Internal
!
  integer                             :: elem
  integer                             :: ii,jj
  integer                             :: propid, propm1
  integer                             :: int_pid
  integer                             :: err_code = 0
  logical, allocatable, dimension(:)  :: knownProp
!
! =================================================================================================
!
! Initialisation
!
  propm1   = 0
!
! =================================================================================================
!
! Calculation
!

  allocate(knownProp(size(fesim%eigenschaften%pbeams,1)))
  knownProp(:) = .FALSE.

  ! Loop over beam2 elements
  do elem=1,size(fesim%elemente%beam2s,1)

    propid = fesim%elemente%beam2s(elem)%pid

    if (propid /= propm1) then

      int_pid = 0
    
      if (fesim%is_pbeam == .true.) then
        do ii=1,size(fesim%eigenschaften%pbeams,1)
          if (propid == fesim%eigenschaften%pbeams(ii)%pid) then

            int_pid  = ii

            if (knownProp(ii) .EQ. .FALSE.) then
              do jj = 1,size(fesim%materialien%mat1s,1)
                if (fesim%eigenschaften%pbeams(ii)%mid == fesim%materialien%mat1s(jj)%mid) then

                  fesim%eigenschaften%pbeams(ii)%lengthWeight = fesim%materialien%mat1s(jj)%rho*fesim%eigenschaften%pbeams(ii)%AA

                  exit
                end if
              end do

              fesim%eigenschaften%pbeams(ii)%weight = 0.d0

              knownProp(ii) = .TRUE.
            end if              

          end if
        end do
      else
        err_code = 1
        goto 9999
      end if

      if (int_pid == 0) then
        write(*,*) 'No property found for beam ', fesim%elemente%beam2s(elem)%eid
        err_code = 2
        goto 9999
      end if 

      propm1 = propid

    end if

    fesim%elemente%beam2s(elem)%int_pid  = int_pid

  end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'beam2_prepareProperties'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine beam2_prepareProperties
