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
!> $Id: sol1_output_hymowi.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_hymowi(Utot)

  use globale_variablen
  use vtk_variablen
  use netz_variablen
  
  implicit none
  
  integer                                            :: err_code=0
  
  integer                                            :: elem

  double precision, dimension(num_dof, num_subcases) :: Utot
  double precision, dimension(num_elements)          :: alphas                  !< Waermeausdehnungskoeffizienten
  character(len = 30)                                :: name
  
  ! Write temperature coefficients for all beam-elements
  alphas = 0.d0
  
  if (is_beam2 == .true.) then
    do elem = 1, size(beam2s,1)
      alphas(beam2s(elem)%eid) = mat1s(pbeams(beam2s(elem)%int_pid)%intMat1ID)%ath
    end do
  
    write(name, '(A11,19x)') 'Alpha_Beam2'
    call vtkoutcelldata(alphas, name)
  end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_hymowi'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_hymowi
