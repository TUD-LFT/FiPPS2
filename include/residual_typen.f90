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
!> $Id: residual_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module residual_typen

  implicit none
  
type residuals_type
  type (initstress_element), allocatable  :: initstrquad8s(:)
  type (initstress_element), allocatable  :: initstr20s(:)
end type residuals_type  
  
! =================================================================================================
!
! Types for storing residual stresses internally
!
type initstress_layer
  double precision, allocatable, dimension(:,:)     :: is_ip
end type initstress_layer

type initstress_element
  type(initstress_layer), allocatable, dimension(:) :: lay
end type initstress_element

contains

  subroutine free_mem_residual(residuals)
  
    implicit none
    
    type(residuals_type) :: residuals
    integer :: ii,jj
    
    if(allocated(residuals%initstrquad8s) .eq. .true.) then
      do ii = 1,size(residuals%initstrquad8s,1)
        do jj = 1,size(residuals%initstrquad8s(ii)%lay,1)
          deallocate(residuals%initstrquad8s(ii)%lay(jj)%is_ip)
        end do
        deallocate(residuals%initstrquad8s(ii)%lay)
      end do
      deallocate(residuals%initstrquad8s)
    end if
    
    if(allocated(residuals%initstr20s) .eq. .true.) then
      do ii = 1,size(residuals%initstr20s,1)
        do jj = 1,size(residuals%initstr20s(ii)%lay,1)
          deallocate(residuals%initstr20s(ii)%lay(jj)%is_ip)
        end do
        deallocate(residuals%initstr20s(ii)%lay)
      end do
      deallocate(residuals%initstr20s)
    end if
  
  end subroutine free_mem_residual

end module residual_typen
