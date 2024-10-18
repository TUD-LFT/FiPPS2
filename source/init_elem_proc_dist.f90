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
!> $Id: init_elem_proc_dist.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine init_elem_proc_dist(fesim,mpisize)

  use fesimulation_typen
  use globale_variablen
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none
  
  type(fe_simulation) :: fesim
  PetscMPIInt :: mpisize
  integer     :: ii, elem, offset
  integer     :: quad8_num, num, rem
  integer     :: lsolid20_num
  integer, dimension(mpisize) :: elem_num_p_proc
  
  if (fesim%is_quad8 .eq. .true.) then
  
    quad8_num = size(fesim%elemente%quad8s,1)
    allocate(mpi_quad8_proc_dist(quad8_num))
  
    mpi_quad8_proc_dist = 0
    
    num = quad8_num / mpisize
    
    rem = MOD(quad8_num,mpisize) 
    
    elem_num_p_proc = num
    
    do ii = 1, rem
    
      elem_num_p_proc(ii) = elem_num_p_proc(ii) + 1
    
    end do
    
    offset = elem_num_p_proc(1)
    do ii = 2, mpisize
      do elem = offset+1, offset+elem_num_p_proc(ii)
        mpi_quad8_proc_dist(elem) = ii-1
      end do
      offset = offset + elem_num_p_proc(ii)
    end do
  
  end if

  if (fesim%is_lsolid20 .eq. .true.) then
  
    lsolid20_num = size(fesim%elemente%lsolid20s,1)
    if (allocated(mpi_lsolid20_proc_dist) .eq. .true.) deallocate(mpi_lsolid20_proc_dist)
    allocate(mpi_lsolid20_proc_dist(lsolid20_num))
  
    mpi_lsolid20_proc_dist = 0
    
    num = lsolid20_num/ mpisize
    
    rem = MOD(lsolid20_num,mpisize) 
    
    elem_num_p_proc = num
    
    do ii = 1, rem
    
      elem_num_p_proc(ii) = elem_num_p_proc(ii) + 1
    
    end do
    
    offset = elem_num_p_proc(1)
    do ii = 2, mpisize
      do elem = offset+1, offset+elem_num_p_proc(ii)
        mpi_lsolid20_proc_dist(elem) = ii-1
      end do
      offset = offset + elem_num_p_proc(ii)
    end do
  
  end if
  
end subroutine init_elem_proc_dist
