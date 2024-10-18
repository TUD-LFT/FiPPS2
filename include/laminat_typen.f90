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
!> $Id: laminat_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module laminat_typen

  implicit none

type laminats_type
  type (lam8_type), allocatable :: lam8s (:)
  type (lam20_type), allocatable :: lam20s(:)
end type laminats_type
! =================================================================================================
!
! Laminate cards
!

!-------------------------------------------------------------------------------
!
! Lam8
!
type lam8_type
  integer               :: lamid            ! Laminate-ID           (Integer>0)
  integer               :: plyid            ! Ply-ID                (Integer>0)
  integer               :: mat8id           ! Mat8-ID               (Integer>0)
  double precision      :: th               ! layer thickness       (Real>0.0)
  double precision      :: angle            ! layer angle           (Real)
end type lam8_type

!-------------------------------------------------------------------------------
!
! Lam20
!
type lam20_type
  integer               :: lamid            ! Laminate-ID           (Integer>0)
  integer               :: plyid            ! Ply-ID                (Integer>0)
  integer,allocatable   :: mat20id(:)       ! Mat20-ID              (Integer>0)
  double precision      :: th               ! layer thickness       (Real>0.0)
  double precision      :: angle            ! layer angle           (Real)
  integer               :: nop              ! number of integration points per layer (Integer)
end type lam20_type

!
! =================================================================================================

contains
  
! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================  
  
  subroutine free_mem_laminate(laminate)
  
    implicit none
    
    type(laminats_type) :: laminate
    integer :: ii
    
    if (allocated(laminate%lam8s)     .eq. .true.) deallocate(laminate%lam8s)
    if (allocated(laminate%lam20s)    .eq. .true.) then
      do ii = 1, size(laminate%lam20s,1)
        deallocate(laminate%lam20s(ii)%mat20id)
      end do
      deallocate(laminate%lam20s)
    end if
  
  end subroutine free_mem_laminate

end module laminat_typen
