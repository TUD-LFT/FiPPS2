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
!> @ details
!> Cell type numbers can be found at "vtkCellType.h"
!> (http://www.vtk.org/doc/release/5.0/html/a03502.html)
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter XX.XX.2009
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter 29.06.2010
!>                                  -> Balkenelemente ergänzt
!
!> $Id: vtkoutelementprop.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
SUBROUTINE vtkoutelementprop(fesim)
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

  INTEGER :: ii
  integer :: offset
  
  integer :: err_code=0
  
  integer, dimension(fesim%num_elements) :: cdata
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

  offset = 0
  IF (ALLOCATED(fesim%elemente%beam2s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%beam2s,1)
      cdata(offset + ii) = fesim%elemente%beam2s(ii)%pid
    END DO
    offset = offset + size(fesim%elemente%beam2s,1)
  END IF
  
  IF (ALLOCATED(fesim%elemente%quad8s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%quad8s,1)
      cdata(offset + ii) = fesim%elemente%quad8s(ii)%pid
    END DO
    offset = offset + size(fesim%elemente%quad8s,1)
  END IF
    
  IF (ALLOCATED(fesim%elemente%lsolid20s) .EQ. .TRUE.) THEN
    DO ii = 1, size(fesim%elemente%lsolid20s,1)
      cdata(offset + ii) = fesim%elemente%lsolid20s(ii)%pid
    END DO
    offset = offset + size(fesim%elemente%lsolid20s,1)
  END IF
  
  CALL vtkoutcelldata_int(cdata, fesim%num_elements, 'Property-ID                   ')
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'vtkoutelementprop'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE vtkoutelementprop
