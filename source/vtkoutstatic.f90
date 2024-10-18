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
SUBROUTINE vtkoutstatic(fesim, Utot, caption, rotations)

  USE fesimulation_typen
  USE vtk_variablen

  IMPLICIT NONE
  TYPE(fe_simulation), INTENT(IN)                        :: fesim
  DOUBLE PRECISION, DIMENSION(fesim%num_dof), INTENT(IN) :: Utot
  CHARACTER(*), INTENT(IN)                               :: caption
  LOGICAL,INTENT(in)                                     :: rotations
  INTEGER                                                :: jj, offset
  

  IF (pointdata .eq. .false.) WRITE(vtkOutFile, '(A11,I8)') 'POINT_DATA ', size(fesim%knoten%nodes,1)

  WRITE(vtkOutFile, '(A)') 'VECTORS ' // TRIM(caption) // ' double'
  
  if (rotations .eq. .false.) then
      offset = 0
  else
      offset = 3
  end if
  
  DO jj = 1, size(fesim%knoten%nodes,1)
      WRITE(vtkOutFile, '(E17.10,X,E17.10,X,E17.10)') Utot((jj-1)*6+1+offset), Utot((jj-1)*6+2+offset), Utot((jj-1)*6+3+offset)
  END DO
 
  celldata = .false.
  pointdata = .true.

END SUBROUTINE vtkoutstatic
