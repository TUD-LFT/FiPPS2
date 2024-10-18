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
subroutine load_add_aeropressure_beam2(fesim, pNode, fac, F_total_aero)

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
  double precision                                  :: length
  double precision, dimension(3)                    :: AB
  type(p2load_type), dimension(:), allocatable      :: origP2loads
  
  fesim%internals%highestLIDtemp = fesim%internals%highestLID
  F_total_struct = 0.d0

  !
  ! Erzeuge die notwendigen P2Loads
  ! 
   
  if (allocated(fesim%lasten%p2loads) .eq. .true.) then
    fesim%internals%indexBeforeAeroP2Load = size(fesim%lasten%p2loads,1)
    allocate(origP2loads(fesim%internals%indexBeforeAeroP2Load))
    
    origP2loads(1:fesim%internals%indexBeforeAeroP2Load) = fesim%lasten%p2loads(1:fesim%internals%indexBeforeAeroP2Load)
    
    deallocate(fesim%lasten%p2loads)
  end if
  
  allocate(fesim%lasten%p2loads(size(fesim%internals%structBeam2IDs,1) + fesim%internals%indexBeforeAeroP2Load))
  
  if (fesim%internals%indexBeforeAeroP2Load .gt. 0) then
    fesim%lasten%p2loads(1:fesim%internals%indexBeforeAeroP2Load) = origP2loads(1:fesim%internals%indexBeforeAeroP2Load)
    deallocate(origP2loads)
  end if
  
  fesim%is_p2load = .true.
  
  do ii = 1, size(fesim%internals%structBeam2IDs,1)
  
    fesim%internals%highestLID = fesim%internals%highestLID + 1
  
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%lid = fesim%internals%highestLID
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%eid1 = fesim%internals%structBeam2IDs(ii)
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%dir = 1
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%thru = .false.
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%eid2 = 0
    
    press = 0.d0
    do jj = 1, 2
      press = press + pNode(fesim%elemente%beam2s(fesim%internals%structBeam2IDs(ii))%nids(jj))
    end do
    
    press = press/2.d0
    
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%pi(1) = press
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%pi(2) = press
    
    AB(:)  = fesim%knoten%nodes(fesim%elemente%beam2s(fesim%internals%structBeam2IDs(ii))%nids(2))%coords(1:3) - fesim%knoten%nodes(fesim%elemente%beam2s(fesim%internals%structBeam2IDs(ii))%nids(1))%coords(1:3)
    
    length = sqrt(dot_product(AB,AB))
    
    F_total_struct = F_total_struct + ABS(press)*length
   
  end do
  
  deallocate(pNode)
  
  do ii = 1, size(fesim%internals%structBeam2IDs,1)
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%pi(1) = fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%pi(1)*F_total_aero(1)/F_total_struct*fac
    fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%pi(2) = fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%pi(2)*F_total_aero(1)/F_total_struct*fac
  end do
  
  
  if (textoutput .eq. .true.) then
    write(*,*) ' F_total_aero  :', F_total_aero(1)
    write(*,*) ' F_total_struct:', F_total_struct
    write(*,*) ' Faktor        :', F_total_aero(1)/F_total_struct
  end if

end subroutine load_add_aeropressure_beam2
