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
!> control subroutine for the creation of internal MPCs
!
!> @details
!> This subrountine controls the convertion of different boundary conditions 
!> to internal MPCs
!> The following conversions are done:
!> - Node-Quad8-Contact -> iMPC 
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 15.10.2012
!
!> $Id: impc_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine impc_control (fesim)

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
type(fe_simulation)         :: fesim
!
! Output
!

!
! inner
!
type (mpc_type),allocatable :: mpcs_temp (:)
type (mpcadd_type),allocatable :: mpcadds_temp(:)

integer                     :: err_code = 0
integer                     :: num_impcs
integer                     :: old_num_mpcs
integer                     :: old_num_mpcadds
integer                     :: ind_offset
integer                     :: ii
!integer                     :: jj

integer, parameter :: numMpcPerContact = 6
integer, parameter :: numMpcPerContactSolidElement = 3

num_impcs = 0

if (fesim%is_contact_node_quad8 .EQ. .TRUE.) then
  do ii = 1,size(fesim%randbedingungen%contact_node_quad8, 1)
    num_impcs = num_impcs + size(fesim%randbedingungen%contact_node_quad8(ii)%nodeIDs, 1) * numMpcPerContact
  end do
end if

if (fesim%is_contact_node_beam2 .EQ. .TRUE.) then
  do ii = 1,size(fesim%randbedingungen%contact_node_beam2, 1)
    num_impcs = num_impcs + size(fesim%randbedingungen%contact_node_beam2(ii)%nodeIDs, 1) * numMpcPerContact
  end do
end if

if (fesim%is_contact_node_lsolid20 .EQ. .TRUE.) then
  do ii = 1,size(fesim%randbedingungen%contact_node_lsolid20, 1)
    num_impcs = num_impcs + size(fesim%randbedingungen%contact_node_lsolid20(ii)%nodeIDs, 1) * numMpcPerContactSolidElement
  end do
end if

if (num_impcs .gt. 0) then

  ind_offset = 1

  if (fesim%is_mpc .eq. .TRUE.) then                                                                    ! Wenn schon MPCs vorhanden sind, muss der Speicher erweitert werden
    old_num_mpcs = size(fesim%randbedingungen%mpcs,1)
    allocate(mpcs_temp(old_num_mpcs))
    mpcs_temp = fesim%randbedingungen%mpcs
    do ii = 1, old_num_mpcs
      deallocate(fesim%randbedingungen%mpcs(ii)%independend)
    end do
    deallocate(fesim%randbedingungen%mpcs)
    allocate(fesim%randbedingungen%mpcs(old_num_mpcs+num_impcs))
    fesim%randbedingungen%mpcs(1:old_num_mpcs) = mpcs_temp
    do ii = 1, old_num_mpcs
      deallocate(mpcs_temp(ii)%independend)
    end do
    deallocate(mpcs_temp)
    ind_offset = old_num_mpcs + 1
    
    old_num_mpcadds = 0
    if (fesim%is_mpcadd .eq. .true.) then
        old_num_mpcadds = size(fesim%randbedingungen%mpcadds,1)
        allocate(mpcadds_temp(old_num_mpcadds))
        mpcadds_temp = fesim%randbedingungen%mpcadds
        deallocate(fesim%randbedingungen%mpcadds)
        allocate(fesim%randbedingungen%mpcadds(old_num_mpcadds+num_impcs))
        fesim%randbedingungen%mpcadds(1:old_num_mpcadds) = mpcadds_temp
    else
        allocate(fesim%randbedingungen%mpcadds(num_impcs))
    end if
    fesim%is_mpcadd = .true.
    do ii = 1,num_impcs
      fesim%randbedingungen%mpcs(old_num_mpcs+ii)%sid = old_num_mpcs+ii
      fesim%randbedingungen%mpcadds(old_num_mpcadds+ii)%scid = -10
      fesim%randbedingungen%mpcadds(old_num_mpcadds+ii)%sid  = old_num_mpcs+ii
    end do
    deallocate(mpcadds_temp)
  else                                                                                          ! Sonst muss der Speicher neu angelegt werden.
    allocate(fesim%randbedingungen%mpcs(num_impcs))
    allocate(fesim%randbedingungen%mpcadds(num_impcs))
    do ii = 1,num_impcs
      fesim%randbedingungen%mpcs(ii)%sid = ii
      fesim%randbedingungen%mpcadds(ii)%scid = -10
      fesim%randbedingungen%mpcadds(ii)%sid  = ii
    end do
    fesim%is_mpc = .true.
    fesim%is_mpcadd = .true.
!    do ii = 1,fesim%num_subcases
!        fesim%lasten%subcases(ii)%mpcaddid = 1
!    end do
  end if

  if (fesim%is_contact_node_quad8 .EQ. .TRUE.) then
    call impc_contact_node_quad8(fesim,ind_offset)
  end if

  if (fesim%is_contact_node_beam2 .EQ. .TRUE.) then
    call impc_contact_node_beam2(fesim,ind_offset)
  end if

  if (fesim%is_contact_node_lsolid20 .EQ. .TRUE.) then
    call impc_contact_node_lsolid20(fesim,ind_offset)
  end if
  
end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'impc_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine impc_control
