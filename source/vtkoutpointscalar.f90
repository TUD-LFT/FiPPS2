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
!> $Id: vtkoutpointscalar.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
SUBROUTINE vtkoutpointscalar(dataP, size, caption)

  USE vtk_variablen

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: size
  INTEGER             :: jj
  CHARACTER(*)        :: caption
  
  double precision, dimension(size), intent(in) :: dataP

  IF (pointdata .eq. .false.) WRITE(vtkOutFile, '(A11,I8)') 'POINT_DATA ', size

  WRITE(vtkOutFile, '(A,A9)') 'SCALARS ' // TRIM(caption), ' double 1'
  WRITE(vtkOutFile, '(A20,A9)') 'LOOKUP_TABLE default'
  
  DO jj = 1, size
    WRITE(vtkOutFile, '(E17.10)') dataP(jj)
  END DO
  
  celldata = .false.
  pointdata = .true.

END SUBROUTINE vtkoutpointscalar