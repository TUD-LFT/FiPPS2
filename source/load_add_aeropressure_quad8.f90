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
subroutine load_add_aeropressure_quad8(fesim, pNode, fac, F_total_aero)

  use fesimulation_typen
  use konstanten

  implicit none
!
! Input+Output
!
  type(fe_simulation)                               :: fesim
  double precision, dimension(:), allocatable       :: pNode
  double precision, intent(in)                      :: fac
  double precision, dimension(:), intent(in)        :: F_total_aero
! 
! Inner
! 
  integer                                           :: ii, jj
  double precision                                  :: press  
  double precision                                  :: F_total_struct
  type(p8load_type), dimension(:), allocatable      :: origP8loads
  
  fesim%internals%highestLIDtemp = fesim%internals%highestLID
  F_total_struct = 0.d0

  !
  ! Erzeuge die notwendigen P8Loads
  ! 
   
  if (allocated(fesim%lasten%p8loads) .eq. .true.) then
    fesim%internals%indexBeforeAeroP8Load = size(fesim%lasten%p8loads,1)
    allocate(origP8loads(fesim%internals%indexBeforeAeroP8Load))
    
    origP8loads(1:fesim%internals%indexBeforeAeroP8Load) = fesim%lasten%p8loads(1:fesim%internals%indexBeforeAeroP8Load)
    
    deallocate(fesim%lasten%p8loads)
  end if
  
  allocate(fesim%lasten%p8loads(size(fesim%internals%structQuad8IDs,1) + fesim%internals%indexBeforeAeroP8Load))
  
  if (fesim%internals%indexBeforeAeroP8Load .gt. 0) then
    fesim%lasten%p8loads(1:fesim%internals%indexBeforeAeroP8Load) = origP8loads(1:fesim%internals%indexBeforeAeroP8Load)
    deallocate(origP8loads)
  end if
  
  fesim%is_p8load = .true.
  
  do ii = 1, size(fesim%internals%structQuad8IDs,1)
  
    fesim%internals%highestLID = fesim%internals%highestLID + 1
  
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%lid = fesim%internals%highestLID
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%eid1 = fesim%internals%structQuad8IDs(ii)
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%cid = 0
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%thru = .false.
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%eid2 = 0
    
    press = 0.d0
    do jj = 1, 8
      press = press + pNode(fesim%elemente%quad8s(fesim%internals%structQuad8IDs(ii))%nids(jj))
    end do
    
    press = -press/8.d0
    
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(1) = press
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(2) = press
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(3) = press
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(4) = press
    
    F_total_struct = F_total_struct + ABS(press)*fesim%elemente%quad8s(fesim%internals%structQuad8IDs(ii))%area
   
  end do
  
  deallocate(pNode)
  
  do ii = 1, size(fesim%internals%structQuad8IDs,1)
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(1) = fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(1)*F_total_aero(1)/F_total_struct*fac
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(2) = fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(2)*F_total_aero(1)/F_total_struct*fac
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(3) = fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(3)*F_total_aero(1)/F_total_struct*fac
    fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(4) = fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%pi(4)*F_total_aero(1)/F_total_struct*fac
  end do
  
  
  if (textoutput .eq. .true.) then
    write(*,*) ' F_total_aero  :', F_total_aero(1)
    write(*,*) ' F_total_struct:', F_total_struct
    write(*,*) ' Faktor        :', F_total_aero(1)/F_total_struct
  end if

end subroutine load_add_aeropressure_quad8
