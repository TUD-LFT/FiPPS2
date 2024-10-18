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
!> VTK-Ergebnisausgabe fuer statische Loesung
!
!> @details
!
!> @author Florian Dexl, TU Dresden, WiMi, 30.09.2021
!
!> $Id: sol1_output_hymowi2.f90 434 2021-01-04 16:50:30Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 434 $
!> $Date: 2021-01-04 17:50:30 +0100 (Mo, 04. Jan 2021) $
!
! =================================================================================================
subroutine sol1_output_vtk(fesim,Utot,scloop)
!
! use
!
  use fesimulation_typen
  use vtk_variablen
  use failure_criteria
  
  implicit none
  
  type(fe_simulation), intent(in)                         :: fesim
  double precision, dimension(fesim%num_dof), intent(in)  :: Utot
  integer,intent(in)                                      :: scloop
  
  integer                                                 :: err_code=0

  integer                                                 :: jj, iLay, maxLay, elem

  double precision, dimension(:,:), allocatable           :: strain, th_strain, mec_strain
  double precision, dimension(:,:), allocatable           :: stress
  double precision, dimension(:,:), allocatable           :: strain_elem, th_strain_elem
  double precision, dimension(:,:), allocatable           :: mec_strain_elem, stress_elem
  double precision, dimension(:,:), allocatable           :: princStrain, princStress
  double precision, dimension(:,:), allocatable           :: princStrain_elem, princStress_elem
  double precision, dimension(:), allocatable             :: nodal_temperatures
  double precision, dimension(:), allocatable             :: beam2temps, quad8temps, lsolid20temps
  double precision, dimension(:), allocatable             :: resFac
  double precision, dimension(:), allocatable             :: tseVec, tseElementalVec
  double precision, dimension(:), allocatable             :: thermal_force
  double precision, dimension(:), allocatable             :: thermal_stress
  integer, dimension(:), allocatable                      :: layNum, failTyp
  character(len = 30)                                     :: name
  double precision, dimension(:), allocatable             :: tstress,tstrain
  double precision                                        :: tResFac
  integer                                                 :: tFailTyp
  double precision, dimension(3)                          :: cent1, cent2, center
  
  ! Routine

  call vtkoutstatic(fesim, Utot, 'Verschiebung', .false.)
  call vtkoutstatic(fesim, Utot, 'Rotation', .true.)

  if (fesim%calculateTSE == .true.) then
    allocate(tseVec(fesim%num_elements))
    tseVec(:) = fesim%ergebnisse%tse
    call vtkoutcelldata(tseVec, fesim%num_elements, 'total_strain_energy')
    deallocate(tseVec)
  end if
  
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
  
    if (fesim%is_beam2 == .true.) then
    
      allocate(stress(fesim%num_elements, 7),thermal_force(fesim%num_elements),thermal_stress(fesim%num_elements))

      stress         = 0.d0
      thermal_force  = 0.d0
      thermal_stress = 0.d0
    
      call beam2_stress(fesim, Utot(:), nodal_temperatures, beam2temps, stress, thermal_force, thermal_stress)

      call vtkoutcelldata(stress(:,1), fesim%num_elements, 'b2_str_dir')
      
      call vtkoutcelldata(stress(:,2), fesim%num_elements, 'b2_str_bndy')
      
      call vtkoutcelldata(stress(:,3), fesim%num_elements, 'b2_str_bndz')
      
      call vtkoutcelldata(stress(:,4), fesim%num_elements, 'b2_str_max')
      
      call vtkoutcelldata(stress(:,5), fesim%num_elements, 'b2_str_min')
      
      allocate(resFac(fesim%num_elements), failTyp(fesim%num_elements), tstress(3), tstrain(3))
      resFac = 0.d0
      
      do elem = 1, size(fesim%elemente%beam2s,1)
      
        tstrain = 0.d0
        tstress = 0.d0
        tstress(1) = stress(fesim%elemente%beam2s(elem)%eid,4)
      
        call getRF_mat1(tstress, &
                      & tstrain, &
                      & fesim%eigenschaften%pbeams(fesim%elemente%beam2s(elem)%int_pid)%intMat1ID, &
                      & resFac(fesim%elemente%beam2s(elem)%eid), &
                      & failTyp(fesim%elemente%beam2s(elem)%eid), &
                      & fesim%versagenskriterien, &
                      & fesim%materialien%mat1s)
                      
        tstress(1) = stress(fesim%elemente%beam2s(elem)%eid,5)
      
        call getRF_mat1(tstress, &
                      & tstrain, &
                      & fesim%eigenschaften%pbeams(fesim%elemente%beam2s(elem)%int_pid)%intMat1ID, &
                      & tResFac, &
                      & tFailTyp, &
                      & fesim%versagenskriterien, &
                      & fesim%materialien%mat1s)
        
        if (tResFac .LT. resFac(fesim%elemente%beam2s(elem)%eid)) then
          resFac(fesim%elemente%beam2s(elem)%eid) = tResFac
          failTyp(fesim%elemente%beam2s(elem)%eid) = tFailTyp
        end if
        
      end do

      call vtkoutcelldata(thermal_force, fesim%num_elements, 'b2_th_force')

      call vtkoutcelldata(thermal_stress, fesim%num_elements, 'b2_th_stress')

      call vtkoutcelldata(resFac, fesim%num_elements, 'b2_ReserveFactor')
      
      call vtkoutcelldata_int(failTyp, size(failTyp,1), 'b2_FailureType')
      
      deallocate(stress, thermal_force, thermal_stress, resFac, failTyp, tstress, tstrain)

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
      
      do elem = 1, size(fesim%elemente%quad8s,1)
      
          cent1(1:3) = 0.5d0*(fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(1))%coords(1:3) + fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(3))%coords(1:3))
          cent2(1:3) = 0.5d0*(fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(2))%coords(1:3) + fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(4))%coords(1:3))
          center(1:3) = 0.5d0*(cent1(1:3) + cent2(1:3))
      
          if (center(1) .LT. fesim%ausgabe%xmin .or. center(1) .GT. fesim%ausgabe%xmax) then
              resFac(fesim%elemente%quad8s(elem)%eid) = 1.d300
              layNum(fesim%elemente%quad8s(elem)%eid) = 0
              failTyp(fesim%elemente%quad8s(elem)%eid) = 0
          end if
      
          if (center(2) .LT. fesim%ausgabe%ymin .or. center(2) .GT. fesim%ausgabe%ymax) then
              resFac(fesim%elemente%quad8s(elem)%eid) = 1.d300
              layNum(fesim%elemente%quad8s(elem)%eid) = 0
              failTyp(fesim%elemente%quad8s(elem)%eid) = 0
          end if
      
          if (center(3) .LT. fesim%ausgabe%zmin .or. center(3) .GT. fesim%ausgabe%zmax) then
              resFac(fesim%elemente%quad8s(elem)%eid) = 1.d300
              layNum(fesim%elemente%quad8s(elem)%eid) = 0
              failTyp(fesim%elemente%quad8s(elem)%eid) = 0
          end if
      
      end do

      if (fesim%calculateElementalTSE .eqv. .true.) then
          allocate(tseElementalVec(fesim%num_elements))
          tseElementalVec    = 0.d0
          call quad8_strainEnergy_control(fesim, Utot(:), nodal_temperatures, quad8temps, tseElementalVec)
          call vtkoutcelldata(tseElementalVec, fesim%num_elements, 'q8_TotalStrainEnergy')
          deallocate(tseElementalVec)
      end if

      call vtkoutcelldata(resFac, fesim%num_elements, 'q8_ReserveFactor')
      
      call vtkoutcelldata_int(layNum, size(layNum,1), 'q8_LayerNumber')
      
      call vtkoutcelldata_int(failTyp, size(failTyp,1), 'q8_FailureType')
      
      deallocate(resFac, layNum, failTyp)
      
#ifndef optitube_output

      call vtkoutcelldata(strain(:,1), fesim%num_elements, 'q8_strain_x')
      
      call vtkoutcelldata(strain(:,2), fesim%num_elements, 'q8_strain_y')
      
      call vtkoutcelldata(strain(:,3), fesim%num_elements, 'q8_strain_xy')
      
      call vtkoutcelldata(strain(:,4), fesim%num_elements, 'q8_strain_1')
      
      call vtkoutcelldata(strain(:,5), fesim%num_elements, 'q8_strain_2')
      
      call vtkoutcelldata(strain(:,6), fesim%num_elements, 'q8_strain_int')
      
      call vtkoutcelldata(stress(:,1), fesim%num_elements, 'q8_stress_x')
      
      call vtkoutcelldata(stress(:,2), fesim%num_elements, 'q8_stress_y')
      
      call vtkoutcelldata(stress(:,3), fesim%num_elements, 'q8_stress_xy')
      
      call vtkoutcelldata(stress(:,4), fesim%num_elements, 'q8_stress_1')
      
      call vtkoutcelldata(stress(:,5), fesim%num_elements, 'q8_stress_2')
      
      call vtkoutcelldata(stress(:,6), fesim%num_elements, 'q8_stress_vMises')

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
            
            write(name, '(A12,I2.2,A4,12x)') 'q8_sig_s_Lay', iLay, '_Top'
            call vtkoutcelldata(stress(:,2), fesim%num_elements, name)
            
            write(name, '(A13,I2.2,A4,12x)') 'q8_sig_ps_Lay', iLay, '_Top'
            call vtkoutcelldata(stress(:,3), fesim%num_elements, name)

            strain  = 0.d0
            stress  = 0.d0
            call quad8_results_elem(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, 2, iLay, .false.)

            write(name, '(A12,I2.2,A4,12x)') 'q8_sig_p_Lay', iLay, '_Bot'
            call vtkoutcelldata(stress(:,1), fesim%num_elements, name)
            
            write(name, '(A12,I2.2,A4,12x)') 'q8_sig_s_Lay', iLay, '_Bot'
            call vtkoutcelldata(stress(:,2), fesim%num_elements, name)
            
            write(name, '(A13,I2.2,A4,12x)') 'q8_sig_ps_Lay', iLay, '_Bot'
            call vtkoutcelldata(stress(:,3), fesim%num_elements, name)
                
        end do
      
      end if
      
#endif
      
!      write(name, '(A16,I2.2,16x)') 'RF_stress_vM_LC_', ii
!      call vtkoutcelldata(stress(:,5), fesim%num_elements, name)
      
      deallocate(stress)
      deallocate(strain)
      
#ifndef optitube_output
      
      ! Knotenergebnisse ausgeben
      
      allocate(strain(fesim%num_nodes, 6))
      allocate(stress(fesim%num_nodes, 6))

      strain = 0.d0
      stress = 0.d0
      
      call quad8_results_node(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, fesim%shellResPos, 0)

      call vtkoutpointscalar(strain(:,1), fesim%num_nodes, 'q8_strain_x')
      
      call vtkoutpointscalar(strain(:,2), fesim%num_nodes, 'q8_strain_y')
      
      call vtkoutpointscalar(strain(:,3), fesim%num_nodes, 'q8_strain_xy')
      
      call vtkoutpointscalar(strain(:,4), fesim%num_nodes, 'q8_strain_1')
      
      call vtkoutpointscalar(strain(:,5), fesim%num_nodes, 'q8_strain_2')
      
      call vtkoutpointscalar(strain(:,6), fesim%num_nodes, 'q8_strain_int')

      call vtkoutpointscalar(stress(:,1), fesim%num_nodes, 'q8_stress_x')
      
      call vtkoutpointscalar(stress(:,2), fesim%num_nodes, 'q8_stress_y')
      
      call vtkoutpointscalar(stress(:,3), fesim%num_nodes, 'q8_stress_xy')
      
      call vtkoutpointscalar(stress(:,4), fesim%num_nodes, 'q8_stress_1')
      
      call vtkoutpointscalar(stress(:,5), fesim%num_nodes, 'q8_stress_2')
      
      call vtkoutpointscalar(stress(:,6), fesim%num_nodes, 'q8_stress_vM')
      
      deallocate(stress, strain)
      
#endif

    end if

    if ( fesim%is_lsolid20  == .true.) then
      
      ! Knotenergebnisse ausgeben

      allocate(th_strain(fesim%num_nodes, 6))
      allocate(strain(fesim%num_nodes, 6))
      allocate(mec_strain(fesim%num_nodes, 6))
      allocate(stress(fesim%num_nodes, 6))
      allocate(th_strain_elem(fesim%num_elements, 6))
      allocate(strain_elem(fesim%num_elements, 6))
      allocate(mec_strain_elem(fesim%num_elements, 6))
      allocate(stress_elem(fesim%num_elements, 6))
      allocate(princStress(fesim%num_nodes,4))
      allocate(princStrain(fesim%num_nodes,4))
      allocate(princStress_elem(fesim%num_elements,4))
      allocate(princStrain_elem(fesim%num_elements,4))
      
      strain = 0.d0
      th_strain = 0.d0
      mec_strain = 0.d0
      stress = 0.d0
      th_strain_elem = 0.d0
      strain_elem = 0.d0
      mec_strain_elem = 0.d0
      stress_elem = 0.d0
      princStress = 0.d0
      princStrain = 0.d0
      princStress_elem = 0.d0
      princStrain_elem = 0.d0

      call lsolid20_strains_stresses(fesim, Utot(:), nodal_temperatures, lsolid20temps, mec_strain, th_strain, strain, stress, princStrain, princStress, mec_strain_elem, th_strain_elem, strain_elem, stress_elem, princStrain_elem, princStress_elem, scloop)

      call vtkoutpointscalar(mec_strain(:,1), fesim%num_nodes, 'l20_me_strain_x')

      call vtkoutpointscalar(mec_strain(:,2), fesim%num_nodes, 'l20_me_strain_y')

      call vtkoutpointscalar(mec_strain(:,3), fesim%num_nodes, 'l20_me_strain_z')

      call vtkoutpointscalar(mec_strain(:,4), fesim%num_nodes, 'l20_me_strain_yz')

      call vtkoutpointscalar(mec_strain(:,5), fesim%num_nodes, 'l20_me_strain_xz')

      call vtkoutpointscalar(mec_strain(:,6), fesim%num_nodes, 'l20_me_strain_xy')

      deallocate(mec_strain)
      
      call vtkoutpointscalar(princStrain(:,1), fesim%num_nodes, 'l20_me_strain_1')

      call vtkoutpointscalar(princStrain(:,2), fesim%num_nodes, 'l20_me_strain_2')

      call vtkoutpointscalar(princStrain(:,3), fesim%num_nodes, 'l20_me_strain_3')

      call vtkoutpointscalar(princStrain(:,4), fesim%num_nodes,  'l20_me_strain_Int')
      
      deallocate(princStrain)
      
      ! if temperature loads are specified
      if ((fesim%is_temperature .EQ. .TRUE.) .OR. (fesim%is_beam2temp .EQ. .TRUE.) .OR. (fesim%is_quad8temp .EQ. .TRUE.) .OR. (fesim%is_lsolid20temp .EQ. .TRUE.)) then
        ! write total strain
        call vtkoutpointscalar(strain(:,1), fesim%num_nodes, 'l20_tt_strain_x')
        
        call vtkoutpointscalar(strain(:,2), fesim%num_nodes, 'l20_tt_strain_y')
        
        call vtkoutpointscalar(strain(:,3), fesim%num_nodes, 'l20_tt_strain_z')
        
        call vtkoutpointscalar(strain(:,4), fesim%num_nodes, 'l20_tt_strain_yz')
        
        call vtkoutpointscalar(strain(:,5), fesim%num_nodes, 'l20_tt_strain_xz')
        
        call vtkoutpointscalar(strain(:,6), fesim%num_nodes, 'l20_tt_strain_xy')

        ! write thermal strain
        call vtkoutpointscalar(th_strain(:,1), fesim%num_nodes, 'l20_th_strain_x')

        call vtkoutpointscalar(th_strain(:,2), fesim%num_nodes, 'l20_th_strain_y')

        call vtkoutpointscalar(th_strain(:,3), fesim%num_nodes, 'l20_th_strain_z')

        call vtkoutpointscalar(th_strain(:,4), fesim%num_nodes, 'l20_th_strain_yz')

        call vtkoutpointscalar(th_strain(:,5), fesim%num_nodes, 'l20_th_strain_xz')

        call vtkoutpointscalar(th_strain(:,6), fesim%num_nodes, 'l20_th_strain_xy')
      end if

      deallocate(strain)
      deallocate(th_strain)
      
      call vtkoutpointscalar(stress(:,1), fesim%num_nodes, 'l20_stress_x')

      call vtkoutpointscalar(stress(:,2), fesim%num_nodes, 'l20_stress_y')

      call vtkoutpointscalar(stress(:,3), fesim%num_nodes, 'l20_stress_z')

      call vtkoutpointscalar(stress(:,4), fesim%num_nodes, 'l20_stress_yz')
      
      call vtkoutpointscalar(stress(:,5), fesim%num_nodes, 'l20_stress_xz')

      call vtkoutpointscalar(stress(:,6), fesim%num_nodes, 'l20_stress_xy')
      
      deallocate(stress)
      
      call vtkoutpointscalar(princStress(:,1), fesim%num_nodes, 'l20_stress_1')

      call vtkoutpointscalar(princStress(:,2), fesim%num_nodes, 'l20_stress_2')

      call vtkoutpointscalar(princStress(:,3), fesim%num_nodes, 'l20_stress_3')

      call vtkoutpointscalar(princStress(:,4), fesim%num_nodes, 'l20_stress_vM')
      
      deallocate(princStress)
      
      call vtkoutcelldata(mec_strain_elem(:,1), fesim%num_elements, 'l20_me_strain_x')
      
      call vtkoutcelldata(mec_strain_elem(:,2), fesim%num_elements, 'l20_me_strain_y')
      
      call vtkoutcelldata(mec_strain_elem(:,3), fesim%num_elements, 'l20_me_strain_z')
      
      call vtkoutcelldata(mec_strain_elem(:,4), fesim%num_elements, 'l20_me_strain_yz')
      
      call vtkoutcelldata(mec_strain_elem(:,5), fesim%num_elements, 'l20_me_strain_xz')
      
      call vtkoutcelldata(mec_strain_elem(:,6), fesim%num_elements, 'l20_me_strain_xy')
      
      deallocate(mec_strain_elem)
      
      call vtkoutcelldata(princStrain_elem(:,1), fesim%num_elements, 'l20_me_strain_1')
      
      call vtkoutcelldata(princStrain_elem(:,2), fesim%num_elements, 'l20_me_strain_2')
      
      call vtkoutcelldata(princStrain_elem(:,3), fesim%num_elements, 'l20_me_strain_3')
      
      call vtkoutcelldata(princStrain_elem(:,4), fesim%num_elements, 'l20_me_strain_Int')
      
      deallocate(princStrain_elem)

      if ((fesim%is_temperature .EQ. .TRUE.) .OR. (fesim%is_beam2temp .EQ. .TRUE.) .OR. (fesim%is_quad8temp .EQ. .TRUE.) .OR. (fesim%is_lsolid20temp .EQ. .TRUE.)) then
        ! write total strain
        call vtkoutcelldata(strain_elem(:,1), fesim%num_elements, 'l20_tt_strain_x')
        
        call vtkoutcelldata(strain_elem(:,2), fesim%num_elements, 'l20_tt_strain_y')
        
        call vtkoutcelldata(strain_elem(:,3), fesim%num_elements, 'l20_tt_strain_z')
        
        call vtkoutcelldata(strain_elem(:,4), fesim%num_elements, 'l20_tt_strain_yz')
        
        call vtkoutcelldata(strain_elem(:,5), fesim%num_elements, 'l20_tt_strain_xz')
        
        call vtkoutcelldata(strain_elem(:,6), fesim%num_elements, 'l20_tt_strain_xy')

        ! write thermal strain
        call vtkoutcelldata(th_strain_elem(:,1), fesim%num_elements, 'l20_th_strain_x')

        call vtkoutcelldata(th_strain_elem(:,2), fesim%num_elements, 'l20_th_strain_y')

        call vtkoutcelldata(th_strain_elem(:,3), fesim%num_elements, 'l20_th_strain_z')

        call vtkoutcelldata(th_strain_elem(:,4), fesim%num_elements, 'l20_th_strain_yz')

        call vtkoutcelldata(th_strain_elem(:,5), fesim%num_elements, 'l20_th_strain_xz')

        call vtkoutcelldata(th_strain_elem(:,6), fesim%num_elements, 'l20_th_strain_xy')
      end if
      
      deallocate(strain_elem)
      deallocate(th_strain_elem)
      
      call vtkoutcelldata(stress_elem(:,1), fesim%num_elements, 'l20_stress_x')

      call vtkoutcelldata(stress_elem(:,2), fesim%num_elements, 'l20_stress_y')

      call vtkoutcelldata(stress_elem(:,3), fesim%num_elements, 'l20_stress_z')

      call vtkoutcelldata(stress_elem(:,4), fesim%num_elements, 'l20_stress_yz')
      
      call vtkoutcelldata(stress_elem(:,5), fesim%num_elements, 'l20_stress_xz')

      call vtkoutcelldata(stress_elem(:,6), fesim%num_elements, 'l20_stress_xy')
      
      deallocate(stress_elem)
      
      call vtkoutcelldata(princStress_elem(:,1), fesim%num_elements, 'l20_stress_1')
      
      call vtkoutcelldata(princStress_elem(:,2), fesim%num_elements, 'l20_stress_2')
      
      call vtkoutcelldata(princStress_elem(:,3), fesim%num_elements, 'l20_stress_3')
      
      call vtkoutcelldata(princStress_elem(:,4), fesim%num_elements, 'l20_stress_vM')
      
      deallocate(princStress_elem)

      allocate(resFac(fesim%num_elements),layNum(fesim%num_elements),failTyp(fesim%num_elements))
      
      resFac = 0.d0
      layNum = -1
      
      call lsolid20_results_elem_lay(fesim, Utot(:), nodal_temperatures, lsolid20temps, resFac, layNum, failTyp, scloop)
      
      call vtkoutcelldata(resFac, fesim%num_elements, 'l20_ReserveFactor')
      
      call vtkoutcelldata_int(layNum, size(layNum,1), 'l20_LayerNumber')
      
      call vtkoutcelldata_int(failTyp, size(failTyp,1), 'l20_FailureType')
      
      deallocate(resFac, layNum, failTyp)

      if (fesim%calculateElementalTSE .eqv. .true.) then
          allocate(tseElementalVec(fesim%num_elements))
          tseElementalVec    = 0.d0
          call lsolid20_strainEnergy_control(fesim, Utot(:), scloop, tseElementalVec)
          call vtkoutcelldata(tseElementalVec, fesim%num_elements, 'l20_TotalStrainEnergy')
          deallocate(tseElementalVec)
      end if
      
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
   write(*,*)                      'sol1_output_vtk'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if  
  
end subroutine sol1_output_vtk
