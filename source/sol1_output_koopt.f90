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
!> $Id: sol1_output_koopt.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_koopt(Utot)

  use globale_variablen
  use netz_variablen
  
  implicit none
  
  integer                                            :: err_code=0
  integer                                            :: outtemp=60
  integer                                            :: ii, elem, kk, pid
  
  double precision, dimension(num_dof, num_subcases) :: Utot
  double precision, dimension(num_elements, 3)       :: epsilon
  double precision, dimension(num_elements, 7)       :: stress
  double precision, dimension(num_nodes, 6)          :: epsilon_node
  double precision, dimension(num_elements, 3)       :: epsilon_node_top
  double precision, dimension(num_elements, 3)       :: epsilon_node_bot
  double precision, dimension(:,:), allocatable      :: stress_dummy
  
  double precision, dimension(:,:,:), allocatable    :: beam_stresses
  double precision, dimension(:,:,:), allocatable    :: shell_strains_max,shell_strains_min
  double precision, dimension(num_subcases)          :: umax
  double precision                                   :: Emat, Gmat, nu
  double precision, dimension(:), allocatable        :: nodal_temperatures
  
  integer                                            :: num_parts, num_pshells
  integer                                            :: pbeamsize
  double precision, dimension(3)                     :: AB
  double precision                                   :: masse
#ifdef vtk_output
  character(len = 30)                                :: name
#endif
  
#ifdef vtk_output
    call vtkoutstatic(Utot, 'Lastfall_', .FALSE.)
#endif  

  allocate(nodal_temperatures(num_nodes))
  nodal_temperatures = 0.d0

  num_parts = 0
  num_pshells = 0
  if (allocated(pshells) .EQ. .TRUE.) then
    num_parts = num_parts + size(pshells,1)
    num_pshells = num_parts
  end if
  if (allocated(pcomps) .EQ. .TRUE.) then
    num_parts = num_parts + size(pcomps,1)
  end if

  pbeamsize = 0
  if (allocated(pbeams)) pbeamsize = size(pbeams,1)

  allocate(beam_stresses(pbeamsize, 4, num_subcases))
  allocate(shell_strains_max(num_parts,3, num_subcases))
  allocate(shell_strains_min(num_parts,3, num_subcases))

  do ii = 1, num_subcases
  
    stress = 0.d0

    if (is_beam2 == .true.) then
    
      call beam2_stress(Utot(:,ii), nodal_temperatures, stress)
    
#ifdef vtk_output
      write(name, '(A12,I2.2,16x)') 'bstr_dir_L1_', ii
      call vtkoutcelldata(stress(:,1), name)
      
      write(name, '(A13,I2.2,16x)') 'bstr_bndy_L1_', ii
      call vtkoutcelldata(stress(:,2), name)
      
      write(name, '(A13,I2.2,16x)') 'bstr_bndz_L1_', ii
      call vtkoutcelldata(stress(:,3), name)
      
      write(name, '(A12,I2.2,16x)') 'bstr_max_L1_', ii
      call vtkoutcelldata(stress(:,4), name)
      
      write(name, '(A12,I2.2,16x)') 'bstr_min_L1_', ii
      call vtkoutcelldata(stress(:,5), name)
#endif

    end if

    if (is_quad8 == .true.) then
    
      allocate(stress_dummy(num_elements, 6))
      stress_dummy = 0.d0

      epsilon_node_top = 0.d0
      !call quad8_results_elem(Utot(:,ii), epsilon_node_top, stress_dummy, 1)
      epsilon_node_bot = 0.d0
      !call quad8_results_elem(Utot(:,ii), epsilon_node_bot, stress_dummy, 2)

#ifdef vtk_output
      epsilon = 0.d0
    
      !call quad8_results_elem(Utot(:,ii), epsilon, stress_dummy, 1)
        
      epsilon_node = 0.d0
      deallocate(stress_dummy)
      allocate(stress_dummy(num_nodes, 6))
      stress_dummy = 0.d0
    
      !call quad8_results_node(Utot(:,ii), epsilon_node, stress_dummy, 0)

      write(name, '(A12,I2.2,16x)') 'strain_x_LCt_', ii
      call vtkoutcelldata(epsilon(:,1), name)
      
      write(name, '(A12,I2.2,16x)') 'strain_y_LCt_', ii
      call vtkoutcelldata(epsilon(:,2), name)
      
      write(name, '(A13,I2.2,16x)') 'strain_xy_LCt_', ii
      call vtkoutcelldata(epsilon(:,3), name)
      
      epsilon = 0.d0
      deallocate(stress_dummy)
      allocate(stress_dummy(num_elements, 6))
      stress_dummy = 0.d0
    
      !call quad8_results_elem(Utot(:,ii), epsilon, stress_dummy, 2)

      write(name, '(A12,I2.2,16x)') 'strain_x_LCb_', ii
      call vtkoutcelldata(epsilon(:,1), name)
      
      write(name, '(A12,I2.2,16x)') 'strain_y_LCb_', ii
      call vtkoutcelldata(epsilon(:,2), name)
      
      write(name, '(A13,I2.2,16x)') 'strain_xy_LCb_', ii
      call vtkoutcelldata(epsilon(:,3), name)

      write(name, '(A12,I2.2,16x)') 'strain_x_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,1), num_nodes, name)
      
      write(name, '(A12,I2.2,16x)') 'strain_y_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,2), num_nodes, name)
      
      write(name, '(A13,I2.2,16x)') 'strain_xy_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,3), num_nodes, name)
#endif
      deallocate(stress_dummy)
    end if
    
    call koopt_sol1_output(stress, epsilon_node_top, epsilon_node_bot, Utot(:,ii), pbeamsize, &
                           umax(ii), beam_stresses(:,:,ii), shell_strains_max(:,:,ii), shell_strains_min(:,:,ii), &
                           num_pshells, num_parts)

  end do  
  
  open(outtemp, file = 'output_FEM.txt', status = 'UNKNOWN')
  
    masse = 0.d0
    do elem = 1,size(beam2s)
      AB = nodes(beam2s(elem)%nids(2))%coords(1:3)-nodes(beam2s(elem)%nids(1))%coords(1:3)
      masse = masse + pbeams(beam2s(elem)%int_pid)%AA*sqrt(dot_product(AB,AB))
    end do
    
    write(outtemp,'(E20.13,20X,A34)') masse, '        Masse                     '
    write(outtemp,*)
  
    do elem = 1,size(pbeams,1)
    
      ! get beam material stiffness
      
      do kk=1,size(mat1s,1)
         if (pbeams(elem)%mid == mat1s(kk)%mid) then
            call mat1_calc_missing_value(kk,Emat,Gmat,nu)
         end if
      end do
    
      do ii = 1, num_subcases
        if (beam_stresses(elem,1,ii) .EQ. -HUGE(1.d0)) beam_stresses(elem,1,ii) = 0.d0
        if (beam_stresses(elem,2,ii) .EQ.  HUGE(1.d0)) beam_stresses(elem,2,ii) = 0.d0
        if (beam_stresses(elem,3,ii) .EQ. -HUGE(1.d0)) beam_stresses(elem,3,ii) = 0.d0
        if (beam_stresses(elem,4,ii) .EQ. -HUGE(1.d0)) beam_stresses(elem,4,ii) = 0.d0
      end do
      write(outtemp,'(I9,A1,I9,A1,A15)') pbeams(elem)%pid, '.', 1, '.', '        Part-ID'
      !write(outtemp,'(2F20.8,A34)') maxval(beam_stresses(elem,1,:)), minval(beam_stresses(elem,2,:)), '        Maximum and minimum stress'
      write(outtemp,'(2E20.13,A34)') maxval(beam_stresses(elem,1,:)), minval(beam_stresses(elem,2,:)), '        Maximum and minimum stress'
      write(outtemp,'(2E20.13,A34)') maxval(beam_stresses(elem,3,:)), maxval(beam_stresses(elem,4,:)), 'max.Axial Force and Bending Moment'
      write(outtemp,'(2F20.8,A35)') maxval(beam_stresses(elem,1,:))/Emat, minval(beam_stresses(elem,2,:))/Emat, '        Maximum and minimum strains'
      write(outtemp,*) 
    end do
    
    do elem = 1,num_parts
    
      if (elem .LT. num_pshells) then
        pid = pshells(elem)%pid
      else
        pid = pcomps(elem)%pid
      end if
    
      write(outtemp,'(I9,A1,I9,A1,A15)') pid, '.', 2, '.', '        Part-ID'
      write(outtemp,'(3F20.8,A23)') minval(shell_strains_min(elem,1,:)), minval(shell_strains_min(elem,2,:)), minval(shell_strains_min(elem,3,:)), '        Minimum strains'
      write(outtemp,'(3F20.8,A23)') maxval(shell_strains_max(elem,1,:)), maxval(shell_strains_max(elem,2,:)), maxval(shell_strains_max(elem,3,:)), '        Maximum strains'
      write(outtemp,*) 
    end do
    
    write(outtemp,'(F20.8,20X,A27)') maxval(umax), '        Maximum deformation'
  
  close(outtemp)
  
  deallocate(beam_stresses)
  
  deallocate(shell_strains_max, shell_strains_min)
  
  deallocate(nodal_temperatures)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_koopt
