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
!> $Id: sol1_output_default.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_default(Utot)

  use globale_variablen
  
  implicit none
  
  integer                                            :: err_code=0
  
  integer                                            :: ii, elem
  
  double precision, dimension(num_dof, num_subcases) :: Utot
  double precision, dimension(num_elements, 3)       :: epsilon
  double precision, dimension(num_elements, 7)       :: stress
  double precision, dimension(num_nodes, 6)          :: epsilon_node
  character(len = 30)                                :: name
  
  call vtkoutstatic(Utot, 'Lastfall_', .FALSE.)

  do ii = 1, num_subcases

    if (is_beam2 == .true.) then
      stress = 0.d0
    
      call beam2_stress(Utot(:,ii), stress)

      write(name, '(A12,I2.2,16x)') 'bstr_dir_LC_', ii
      call vtkoutcelldata(stress(:,1), name)
      
      write(name, '(A13,I2.2,16x)') 'bstr_bndy_LC_', ii
      call vtkoutcelldata(stress(:,2), name)
      
      write(name, '(A13,I2.2,16x)') 'bstr_bndz_LC_', ii
      call vtkoutcelldata(stress(:,3), name)
      
      write(name, '(A12,I2.2,16x)') 'bstr_max_LC_', ii
      call vtkoutcelldata(stress(:,4), name)
      
      write(name, '(A12,I2.2,16x)') 'bstr_min_LC_', ii
      call vtkoutcelldata(stress(:,5), name)

    end if

    if (is_quad8 == .true.) then

      epsilon = 0.d0
    
      call quad8_strains(Utot(:,ii), epsilon, 1)
      
      epsilon_node = 0.d0
    
      call quad8_strains_node(Utot(:,ii), epsilon_node, 0)

      write(name, '(A12,I2.2,16x)') 'strain_x_LC_', ii
      call vtkoutcelldata(epsilon(:,1), name)
      
      write(name, '(A12,I2.2,16x)') 'strain_y_LC_', ii
      call vtkoutcelldata(epsilon(:,2), name)
      
      write(name, '(A13,I2.2,16x)') 'strain_xy_LC_', ii
      call vtkoutcelldata(epsilon(:,3), name)

      write(name, '(A12,I2.2,16x)') 'strain_x_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,1), num_nodes, name)
      
      write(name, '(A12,I2.2,16x)') 'strain_y_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,2), num_nodes, name)
      
      write(name, '(A13,I2.2,16x)') 'strain_xy_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,3), num_nodes, name)
      
      write(name, '(A12,I2.2,16x)') 'strain_1_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,4), num_nodes, name)
      
      write(name, '(A12,I2.2,16x)') 'strain_2_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,5), num_nodes, name)
      
      write(name, '(A12,I2.2,16x)') 'strain_int_LC_', ii
      call vtkoutpointscalar(epsilon_node(:,6), num_nodes, name)

    end if

  end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_default'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_default
