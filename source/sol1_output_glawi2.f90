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
!
!> $Id: sol1_output_glawi2.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
subroutine sol1_output_glawi2(fesim,Utot,scloop,cdi,aeroConverged)
! =================================================================================================
!
!       Header:         
!
!       Content:        
!
!       Input:          
!
!       Output:         
!
!       Calls:          
!
!       Called by:      
!
!       Author:         
!
! =================================================================================================
!
! use
!
  use fesimulation_typen
  use vtk_variablen
  use failure_criteria
  
  implicit none
  
  type(fe_simulation)                                :: fesim
  double precision, dimension(fesim%num_dof)         :: Utot
  integer,intent(in)                                 :: scloop
  double precision, dimension(:)                     :: cdi
  logical                                            :: aeroConverged
  
  integer                                            :: err_code=0
  
  integer                                            :: ii, elem
  
  double precision, dimension(:,:), allocatable      :: strain
  double precision, dimension(:,:), allocatable      :: stress
  double precision, dimension(:), allocatable        :: nodal_temperatures
  double precision, dimension(:), allocatable        :: beam2temps, quad8temps, lsolid20temps
  double precision, dimension(:), allocatable        :: resFac
  integer, dimension(:), allocatable                 :: layNum, failTyp
  integer, parameter                                 :: outfem = 55
  double precision                                   :: masse
  
  if (scloop .eq. 1) then
    OPEN(outfem, file = 'output_FEM.txt', STATUS='UNKNOWN')
    
    masse = 0.d0

    do ii = 1, size(fesim%eigenschaften%pcomps,1)
        masse = masse + fesim%eigenschaften%pcomps(ii)%weight
    end do

    do ii = 1, size(fesim%eigenschaften%pshells,1)
        masse = masse + fesim%eigenschaften%pshells(ii)%weight
    end do

    write(outfem,'(A21,E25.18)') 'Ges.-masse:          ', masse
    write(outfem,'(A21,I25)')    'max. Lastfallanzahl: ', fesim%num_subcases
  else
    OPEN(outfem, file = 'output_FEM.txt', STATUS='OLD', POSITION='APPEND')
  end if
  
  if (aeroConverged .eq. .true.) then
    allocate(nodal_temperatures(fesim%num_nodes))
    nodal_temperatures = 0.d0
  
    if (fesim%is_temperature == .true.) then
      call get_node_temperatures(fesim,fesim%lasten%subcases(scloop)%loadid, nodal_temperatures)
    end if

    allocate(beam2temps(size(fesim%elemente%beam2s,1)), quad8temps(size(fesim%elemente%quad8s,1)), lsolid20temps(size(fesim%elemente%lsolid20s,1)))
    beam2temps = 0.d0
    quad8temps = 0.d0
    lsolid20temps = 0.d0

    ! Get elemental temperatures at beam2 elements
    if (fesim%is_beam2temp .eqv. .true.) call get_beam2_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, beam2temps)
    ! Get elemental temperatures at quad8 elements
    if (fesim%is_quad8temp .eqv. .true.) call get_quad8_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, quad8temps)
    ! Get elemental temperatures at lsolid20 elements
    if (fesim%is_lsolid20temp .eqv. .true.) call get_lsolid20_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, lsolid20temps)

    if (fesim%is_quad8 == .true.) then
      
      
      allocate(resFac(fesim%num_elements),layNum(fesim%num_elements),failTyp(fesim%num_elements))
      
      resFac = 0.d0
      layNum = -1
      
      call quad8_results_elem_lay(fesim, Utot(:), nodal_temperatures, quad8temps, resFac, layNum, failTyp)
      
      ! Elementergebnisse ausgeben
      
      allocate(strain(fesim%num_elements, 6))
      allocate(stress(fesim%num_elements, 6))

      strain  = 0.d0
      stress  = 0.d0
    
      call quad8_results_elem(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, fesim%shellResPos, 0, .true.)
      
      do elem = 1, size(fesim%elemente%quad8s,1)
      
        if (fesim%elemente%quad8s(elem)%propType == 1) then
        
          call getRF_mat1(stress(fesim%elemente%quad8s(elem)%eid,1:3), &
                & strain(fesim%elemente%quad8s(elem)%eid,1:3), &
                & fesim%eigenschaften%pshells(fesim%elemente%quad8s(elem)%int_pid)%intMat1ID, &
                & resFac(fesim%elemente%quad8s(elem)%eid), &
                & failTyp(fesim%elemente%quad8s(elem)%eid), &
                & fesim%versagenskriterien, &
                & fesim%materialien%mat1s)
          layNum(fesim%elemente%quad8s(elem)%eid) = 0
                  
        end if
      
      end do
      
      write(outfem,'(A21,E25.18)') 'min. Reservefactor:  ', minval(resFac(fesim%elemente%quad8s(1)%eid:fesim%elemente%quad8s(size(fesim%elemente%quad8s,1))%eid))
      write(outfem,'(A21,E25.18)') 'ind. Drag-Coeff:     ', cdi(1)
      write(outfem,'(A21,E25.18)') 'max. princ. stress:  ', maxval(stress(fesim%elemente%quad8s(1)%eid:fesim%elemente%quad8s(size(fesim%elemente%quad8s,1))%eid,4))
      
      deallocate(stress, strain)
      deallocate(resFac, layNum, failTyp)

    end if
  
    deallocate(nodal_temperatures, beam2temps, quad8temps, lsolid20temps)
  end if
  
  close(outfem)
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_glawi2'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_glawi2
