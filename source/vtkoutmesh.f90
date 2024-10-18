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
!> This Subroutine creates the mesh for the ParaView output
!
!> @details
!> Cell type numbers can be found at "vtkCellType.h"
!> (http://www.vtk.org/doc/release/5.0/html/a03502.html)
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter XX.XX.2009
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter 29.06.2010
!>                                      -> Balkenelemente ergänzt 
!
!> $Id: vtkoutmesh.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
SUBROUTINE vtkoutmesh(fesim)
!
! Use
!
  USE vtk_variablen
  USE fesimulation_typen
!
! =================================================================================================
!
  IMPLICIT NONE
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Data types
!
  type(fe_simulation) :: fesim

  INTEGER :: ii,jj
  INTEGER :: cell_number, cell_data_number
  
  integer :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
! 
  
  ! Knoten
  
  WRITE(vtkOutFile, '(A7,I8,A7)') 'POINTS ', size(fesim%knoten%nodes,1), ' double'
  
  DO ii = 1, size(fesim%knoten%nodes,1)
      WRITE(vtkOutFile, '(E24.17,X,E24.17,X,E24.17)') fesim%knoten%nodes(ii)%coords(1), fesim%knoten%nodes(ii)%coords(2), fesim%knoten%nodes(ii)%coords(3)
  END DO
  
  ! Anzahl Zellen und Zelldaten
  
  cell_number      = 0
  cell_data_number = 0
  
  IF (ALLOCATED(fesim%elemente%beam2s) .EQ. .TRUE.) THEN 
    cell_number      = cell_number + size(fesim%elemente%beam2s,1)
    cell_data_number = cell_data_number + size(fesim%elemente%beam2s,1) * 3
  END IF
  
  IF (ALLOCATED(fesim%elemente%quad8s) .EQ. .TRUE.) THEN 
    cell_number      = cell_number + size(fesim%elemente%quad8s,1)
    cell_data_number = cell_data_number + size(fesim%elemente%quad8s,1) * 9
  END IF
  
  IF (ALLOCATED(fesim%elemente%lsolid20s) .EQ. .TRUE.) THEN
    cell_number      = cell_number + size(fesim%elemente%lsolid20s,1)
    cell_data_number = cell_data_number + size(fesim%elemente%lsolid20s,1) * 21
  END IF
  
  WRITE(vtkOutFile, '(A6,I8,X,I8)') 'CELLS ', cell_number, cell_data_number
  
  ! Zellen schreiben
  
  IF (ALLOCATED(fesim%elemente%beam2s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%beam2s,1)
      WRITE(vtkOutFile, '(2(I10,X),I10)') 2, fesim%elemente%beam2s(ii)%nids(1)-1, fesim%elemente%beam2s(ii)%nids(2)-1
    END DO
  END IF
  
  IF (ALLOCATED(fesim%elemente%quad8s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%quad8s,1)
      WRITE(vtkOutFile, '(8(I10,X),I10)') 8, fesim%elemente%quad8s(ii)%nids(1)-1, fesim%elemente%quad8s(ii)%nids(2)-1, fesim%elemente%quad8s(ii)%nids(3)-1, fesim%elemente%quad8s(ii)%nids(4)-1, fesim%elemente%quad8s(ii)%nids(5)-1, fesim%elemente%quad8s(ii)%nids(6)-1, fesim%elemente%quad8s(ii)%nids(7)-1, fesim%elemente%quad8s(ii)%nids(8)-1
    END DO
  END IF

  IF (ALLOCATED(fesim%elemente%lsolid20s) .EQ. .TRUE.) THEN
    DO ii = 1, size(fesim%elemente%lsolid20s,1)
      WRITE(vtkOutFile, '(20(I10,X),I10)') 20, (fesim%elemente%lsolid20s(ii)%nids(jj)-1,jj=1,20)
    END DO
  END IF
 
  WRITE(vtkOutFile, '(A11,I8)') 'CELL_TYPES ', cell_number
  
  ! Zelltypen zuordnen
  
  IF (ALLOCATED(fesim%elemente%beam2s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%beam2s,1)
      WRITE(vtkOutFile, '(I10)') 3
    END DO
  END IF
  
  IF (ALLOCATED(fesim%elemente%quad8s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%quad8s,1)
      WRITE(vtkOutFile, '(I10)') 23
    END DO
  END IF
  
  IF (ALLOCATED(fesim%elemente%lsolid20s) .EQ. .TRUE.) THEN
    DO ii = 1, size(fesim%elemente%lsolid20s,1)
      ! Cell type 25: quadratic hexahedron
      WRITE(vtkOutFile, '(I10)') 25
    END DO
  END IF
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'vtkoutmesh'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE vtkoutmesh
