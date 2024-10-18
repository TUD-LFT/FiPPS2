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
!> $Id: vtkouteig.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
SUBROUTINE vtkouteig(fesim, eigenvalnum, eigenvalue, eigenvector)
!
! Use
!
  use fesimulation_typen
  USE vtk_variablen
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
!
! Input
!
  type(fe_simulation)                                   :: fesim
  integer, intent(in)                                   :: eigenvalnum
  double precision, intent(in)                          :: eigenvalue
  double precision, dimension(fesim%num_dof), intent(in):: eigenvector
!
! Internal
!
  integer                                               :: jj
  
  integer                                               :: err_code=0
  
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
  IF (pointdata .eq. .false.) WRITE(vtkOutFile, '(A11,I8)') 'POINT_DATA ', size(fesim%knoten%nodes,1)
  
! Eigenform herausschreiben
  WRITE(vtkOutFile, '(A20,I2.2,A7)') 'VECTORS Eigenvektor_', eigenvalnum,' double'
  DO jj = 1, size(fesim%knoten%nodes,1)
      WRITE(vtkOutFile, '(E17.10,X,E17.10,X,E17.10)') eigenvector((jj-1)*6+1), eigenvector((jj-1)*6+2), eigenvector((jj-1)*6+3)
  END DO

! Eigenwert herausschreiben
  WRITE(vtkOutFile, '(A10,I8)') 'CELL_DATA ', fesim%num_elements
  WRITE(vtkOutFile, '(A19,I2.2,A9)') 'SCALARS Eigenvalue_', eigenvalnum,' double 1'
  WRITE(vtkOutFile, '(A20)')         'LOOKUP_TABLE default'
  DO jj = 1, fesim%num_elements
      WRITE(vtkOutFile, '(E18.10E3)') eigenvalue
  END DO
  
  celldata = .true.
  pointdata = .false.

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'vtkouteig'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE vtkouteig
