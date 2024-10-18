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
!> $Id: vtkoutcelldata.f90 362 2018-08-07 09:18:14Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 362 $
!> $Date: 2018-08-07 11:18:14 +0200 (Di, 07. Aug 2018) $
!
! =================================================================================================
SUBROUTINE vtkoutcelldata_vector(cdata, num, caption)

  USE vtk_variablen

  IMPLICIT NONE
  INTEGER :: jj
  CHARACTER(*) :: caption
  
  integer :: num
  double precision, dimension(num,3), intent(in) :: cdata

  IF (celldata .EQ. .FALSE.) WRITE(vtkOutFile, '(A10,I8)') 'CELL_DATA ', num

  WRITE(vtkOutFile, '(A,A9)') 'VECTORS ' // TRIM(caption), ' double'
  
  DO jj = 1, num
      WRITE(vtkOutFile, '(E18.10E3,X,E18.10E3,X,E18.10E3)') cdata(jj,1), cdata(jj,2), cdata(jj,3)
  END DO
  
  celldata = .true.
  pointdata = .false.

END SUBROUTINE vtkoutcelldata_vector
