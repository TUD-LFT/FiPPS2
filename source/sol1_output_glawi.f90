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
!> $Id: sol1_output_glawi.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_glawi(fesim,Utot,scloop)
!
! use
!
  use fesimulation_typen
  use vtk_variablen
  use failure_criteria
  use konstanten
  
  implicit none
  
  type(fe_simulation)                                :: fesim
  integer,intent(in)                                 :: scloop
  
  integer                                            :: err_code=0
  
  integer                                            :: jj, iLay, maxLay, elem
  
  double precision, dimension(fesim%num_dof)         :: Utot
  double precision, dimension(:,:), allocatable      :: strain
  double precision, dimension(:,:), allocatable      :: stress
  double precision, dimension(:), allocatable        :: nodal_temperatures
  double precision, dimension(:), allocatable        :: beam2temps, quad8temps, lsolid20temps
  double precision, dimension(:), allocatable        :: resFac
  double precision, dimension(:), allocatable        :: thermal_force
  double precision, dimension(:), allocatable        :: thermal_stress
  integer, dimension(:), allocatable                 :: layNum, failTyp
  double precision, dimension(:), allocatable        :: angle                   !< Richtung der Hauptdehnung (Winkel zwischen dem Elementkoordinatensystem und der 1. Hauptdehnung
  character(len = 30)                                :: name
  
  call vtkoutstatic(fesim, Utot, 'Verschiebung', .false.)
  
  allocate(nodal_temperatures(fesim%num_nodes))
  nodal_temperatures = 0.d0
  
  allocate(beam2temps(size(fesim%elemente%beam2s,1)))
  allocate(quad8temps(size(fesim%elemente%quad8s,1)))
  allocate(lsolid20temps(size(fesim%elemente%lsolid20s,1)))
  beam2temps = 0.d0
  quad8temps = 0.d0
  lsolid20temps = 0.d0
  
    if (fesim%is_temperature == .true.) then
      call get_node_temperatures(fesim,fesim%lasten%subcases(scloop)%loadid, nodal_temperatures)
    end if

    ! Get elemental temperatures at beam2 elements
    if (fesim%is_beam2temp .eqv. .true.) call get_beam2_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, beam2temps)
    ! Get elemental temperatures at quad8 elements
    if (fesim%is_quad8temp .eqv. .true.) call get_quad8_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, quad8temps)
    ! Get elemental temperatures at lsolid20 elements
    if (fesim%is_lsolid20temp .eqv. .true.) call get_lsolid20_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, lsolid20temps)
  
    if (fesim%is_beam2 == .true.) then
    
      allocate(stress(fesim%num_elements, 7),thermal_force(fesim%num_elements),thermal_stress(fesim%num_elements))
      
      stress = 0.d0
      thermal_force  = 0.d0
      thermal_stress = 0.d0
    
      call beam2_stress(fesim, Utot(:), nodal_temperatures, beam2temps, stress, thermal_force, thermal_stress)

      call vtkoutcelldata(stress(:,1), fesim%num_elements, 'b2_str_dir')
      
      deallocate(stress,thermal_force,thermal_stress)

    end if

    if (fesim%is_quad8 == .true.) then
    
      ! Elementergebnisse ausgeben
      
      allocate(strain(fesim%num_elements, 6))
      allocate(stress(fesim%num_elements, 6))

      strain  = 0.d0
      stress  = 0.d0
    
      call quad8_results_elem(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, fesim%shellResPos, 0, .false.)
      
      
      allocate(resFac(fesim%num_elements),layNum(fesim%num_elements),failTyp(fesim%num_elements))
      
      resFac = 0.d0
      layNum = -1
      
      call quad8_results_elem_lay(fesim, Utot(:), nodal_temperatures, quad8temps, resFac, layNum, failTyp)
      
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
      
      call vtkoutcelldata(resFac, fesim%num_elements, 'q8_ReserveFactor')
      
      call vtkoutcelldata_int(layNum, size(layNum,1), 'q8_LayerNumber')
      
      call vtkoutcelldata_int(failTyp, size(failTyp,1), 'q8_FailureType')
      
      deallocate(resFac, layNum, failTyp)
      
      call vtkoutcelldata(strain(:,4), fesim%num_elements, 'q8_strain_1')
      
      call vtkoutcelldata(strain(:,5), fesim%num_elements, 'q8_strain_2')
      
      call vtkoutcelldata(stress(:,4), fesim%num_elements, 'q8_stress_1')
      
      call vtkoutcelldata(stress(:,5), fesim%num_elements, 'q8_stress_2')
      
      ! Hauptdehnungsrichtung
      ! Für das GLARE sind die Fatiguewerte für Winkel kleiner 20° vom Winkel nahezu unabhängig.
      
      allocate(angle(fesim%num_elements))
      
      angle = 0.d0
      
      do elem = 1, size(fesim%elemente%quad8s,1)
      
        if (fesim%elemente%quad8s(elem)%propType == 2) then
      
          angle(fesim%elemente%quad8s(elem)%eid) = 0.5d0 * atan(strain(fesim%elemente%quad8s(elem)%eid,3) / (strain(fesim%elemente%quad8s(elem)%eid,1) - strain(fesim%elemente%quad8s(elem)%eid,2)))
        
          angle(fesim%elemente%quad8s(elem)%eid) = angle(fesim%elemente%quad8s(elem)%eid) * 180.d0 / PI
        
        end if
      
      end do
      
      write(name, '(A14,16x)') 'q8_pstrain_ang'
      call vtkoutcelldata(angle, fesim%num_elements, name)
      
      deallocate(angle)

      if (fesim%is_pcomp .eq. .true.) then
      
        maxLay = 0
        do jj = 1, size(fesim%eigenschaften%pcomps,1) 
          if (maxLay .lt. fesim%eigenschaften%pcomps(jj)%lay) then
            maxLay = fesim%eigenschaften%pcomps(jj)%lay
          end if
        end do
        
        do iLay = 1, maxLay
            strain  = 0.d0
            stress  = 0.d0
            call quad8_results_elem(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, 1, iLay, .false.)

            write(name, '(A12,I2.2,A4,12x)') 'q8_sig_p_Lay', iLay, '_Top'
            call vtkoutcelldata(stress(:,1), fesim%num_elements, name)

            strain  = 0.d0
            stress  = 0.d0
            call quad8_results_elem(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, 2, iLay, .false.)

            write(name, '(A12,I2.2,A4,12x)') 'q8_sig_p_Lay', iLay, '_Bot'
            call vtkoutcelldata(stress(:,1), fesim%num_elements, name)
                
        end do
      
      end if
      
      deallocate(stress)
      deallocate(strain)

    end if
  
  deallocate(nodal_temperatures, beam2temps, quad8temps, lsolid20temps)
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_glawi'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_glawi
