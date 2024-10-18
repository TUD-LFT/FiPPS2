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
! =================================================================================================
SUBROUTINE vtkoutpointvektor(pdata, size, caption)

  USE fesimulation_typen
  USE vtk_variablen

  IMPLICIT NONE
  INTEGER, INTENT(IN)                                    :: size
  DOUBLE PRECISION, DIMENSION(size,3), INTENT(IN)        :: pdata
  CHARACTER(*), INTENT(IN)                               :: caption
  INTEGER                                                :: jj
  

  IF (pointdata .eq. .false.) WRITE(vtkOutFile, '(A11,I8)') 'POINT_DATA ', size

  WRITE(vtkOutFile, '(A)') 'VECTORS ' // TRIM(caption) // ' double'
  
  DO jj = 1, size
      WRITE(vtkOutFile, '(E17.10,X,E17.10,X,E17.10)') pdata(jj,1), pdata(jj,2), pdata(jj,3)
  END DO
 
  celldata = .false.
  pointdata = .true.

END SUBROUTINE vtkoutpointvektor
