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
subroutine load_add_aeropressure(fesim, scloop, pAeroElem, F_total_aero)

  use fesimulation_typen
  use konstanten

  implicit none
!
! Input+Output
!
  type(fe_simulation)                               :: fesim
  integer, intent(in)                               :: scloop
  double precision, dimension(:,:), intent(in)      :: pAeroElem
  double precision, dimension(:), intent(in)        :: F_total_aero
! 
! Inner
! 
  integer                                           :: ii
  integer                                           :: loadID
  double precision                                  :: fac
  double precision, dimension(:), allocatable       :: pNode
  type(load_type), dimension(:), allocatable        :: origLoads
  
  interface
    subroutine load_add_aeropressure_beam2(fesim, pNode, fac, F_total_aero)
        use fesimulation_typen

        type(fe_simulation)                               :: fesim
        double precision, dimension(:), allocatable       :: pNode
        double precision, intent(in)                      :: fac
        double precision, dimension(:), intent(in)        :: F_total_aero
    end subroutine load_add_aeropressure_beam2

    subroutine load_add_aeropressure_quad8(fesim, pNode, fac, F_total_aero)
        use fesimulation_typen

        type(fe_simulation)                               :: fesim
        double precision, dimension(:), allocatable       :: pNode
        double precision, intent(in)                      :: fac
        double precision, dimension(:), intent(in)        :: F_total_aero
    end subroutine load_add_aeropressure_quad8
  end interface

  !
  ! Bestimmen der Drücke an jedem Strukturknoten
  !

  loadID = fesim%lasten%subcases(scloop)%loadid
    
  allocate(pNode(fesim%num_nodes))
  
  pNode = 0.d0
  
  do ii = 1, size(fesim%internals%aeroElem2structNode,1)
  
    if (pNode(fesim%internals%aeroElem2structNode(ii)%nodeID) .ne. 0.d0) then
      write(*,*) 'Doppelte Zuweisung auf einen Strukturknoten'
    end if
  
    pNode(fesim%internals%aeroElem2structNode(ii)%nodeID) = pAeroElem(fesim%internals%aeroElem2structNode(ii)%elemID,1)
  
  end do

  fac = fesim%lasten%subcases(scloop)%aeroloadFac

  if (fesim%is_aeroload3d .eq. .true.) then
    call load_add_aeropressure_quad8(fesim, pNode, fac, F_total_aero)
  end if
  
  if (fesim%is_aeroload2d .eq. .true.) then
    call load_add_aeropressure_beam2(fesim, pNode, fac, F_total_aero)
  end if
  
  !
  ! Erzeuge die notwendigen Loads
  ! 
  
  if (allocated(fesim%lasten%loads) .eq. .true.) then
    fesim%internals%indexBeforeAeroLoad = size(fesim%lasten%loads,1)
    allocate(origLoads(fesim%internals%indexBeforeAeroLoad))
    
    origLoads(1:fesim%internals%indexBeforeAeroLoad) = fesim%lasten%loads(1:fesim%internals%indexBeforeAeroLoad)
    
    deallocate(fesim%lasten%loads)
  end if
  
  allocate(fesim%lasten%loads(size(fesim%internals%structBeam2IDs,1) + size(fesim%internals%structQuad8IDs,1) + fesim%internals%indexBeforeAeroLoad))
  
  if (fesim%internals%indexBeforeAeroLoad .gt. 0) then
    fesim%lasten%loads(1:fesim%internals%indexBeforeAeroLoad) = origloads(1:fesim%internals%indexBeforeAeroLoad)
    deallocate(origloads)
  end if
  
  fesim%is_load = .true.
  
  do ii = 1, size(fesim%internals%structBeam2IDs,1)
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%lidi = fesim%lasten%p2loads(fesim%internals%indexBeforeAeroP2Load+ii)%lid
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%lcid = loadID
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%sfaci = 1.d0
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%sfac = 1.d0
  end do
  
  do ii = 1, size(fesim%internals%structQuad8IDs,1)
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%lidi = fesim%lasten%p8loads(fesim%internals%indexBeforeAeroP8Load+ii)%lid
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%lcid = loadID
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%sfaci = 1.d0
    fesim%lasten%loads(fesim%internals%indexBeforeAeroLoad+ii)%sfac = 1.d0
  end do

end subroutine load_add_aeropressure
