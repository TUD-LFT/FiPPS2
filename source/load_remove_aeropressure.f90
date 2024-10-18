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
subroutine load_remove_aeropressure(fesim)

  use fesimulation_typen

  implicit none
!
! Input+Output
!
  type(fe_simulation)                               :: fesim
! 
! Inner
! 
  type(p2load_type), dimension(:), allocatable      :: origP2loads
  type(p8load_type), dimension(:), allocatable      :: origP8loads
  type(load_type), dimension(:), allocatable        :: origLoads
  
  if (fesim%internals%indexBeforeAeroP2Load .eq. 0) then                        ! vorher gab es keine p2loads
    fesim%is_p2load = .false.
    if (allocated(fesim%lasten%p2loads)) deallocate(fesim%lasten%p2loads)
  else                                                                          ! sonst
    
    allocate(origP2loads(fesim%internals%indexBeforeAeroP2Load))
    origP2loads(1:fesim%internals%indexBeforeAeroP2Load) = fesim%lasten%p2loads(1:fesim%internals%indexBeforeAeroP2Load)
    deallocate(fesim%lasten%p2loads)
    
    allocate(fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load))
    fesim%lasten%p2loads(1:fesim%internals%indexBeforeAeroP2Load) = origP2loads(1:fesim%internals%indexBeforeAeroP2Load)
    deallocate(origP2loads)
    
    fesim%internals%indexBeforeAeroP2Load = 0
    fesim%is_p2load = .true.
  
  end if
  
  if (fesim%internals%indexBeforeAeroP8Load .eq. 0) then                        ! vorher gab es keine p8loads
    fesim%is_p8load = .false.
    if (allocated(fesim%lasten%p8loads)) deallocate(fesim%lasten%p8loads)
  else                                                                          ! sonst
    
    allocate(origP8loads(fesim%internals%indexBeforeAeroP8Load))
    origP8loads(1:fesim%internals%indexBeforeAeroP8Load) = fesim%lasten%p8loads(1:fesim%internals%indexBeforeAeroP8Load)
    deallocate(fesim%lasten%p8loads)
    
    allocate(fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load))
    fesim%lasten%p8loads(1:fesim%internals%indexBeforeAeroP8Load) = origP8loads(1:fesim%internals%indexBeforeAeroP8Load)
    deallocate(origP8loads)
    
    fesim%internals%indexBeforeAeroP8Load = 0
    fesim%is_p8load = .true.
  
  end if
  
  if (fesim%internals%indexBeforeAeroLoad .eq. 0) then                          ! vorher gab es keine loads
    fesim%is_load = .false.
    if (allocated(fesim%lasten%loads)) deallocate(fesim%lasten%loads)
  else                                                                          ! sonst
    
    allocate(origLoads(fesim%internals%indexBeforeAeroLoad))
    origLoads(1:fesim%internals%indexBeforeAeroLoad) = fesim%lasten%loads(1:fesim%internals%indexBeforeAeroLoad)
    deallocate(fesim%lasten%loads)
    
    allocate(fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad))
    fesim%lasten%loads(1:fesim%internals%indexBeforeAeroLoad) = origLoads(1:fesim%internals%indexBeforeAeroLoad)
    deallocate(origLoads)
    
    fesim%internals%indexBeforeAeroLoad = 0
    fesim%is_load = .true.
  
  end if
  
  fesim%internals%highestLID = fesim%internals%highestLIDtemp

end subroutine load_remove_aeropressure
