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
!> $Id: sol2_output.f90 414 2019-10-17 11:45:47Z DOM+ahauffe $
!> $Author: DOM+ahauffe $
!> $Revision: 414 $
!> $Date: 2019-10-17 13:45:47 +0200 (Do, 17. Okt 2019) $
!
! =================================================================================================
SUBROUTINE sol2_output_control(fesim, sc, eigenvalues, eigenvectors, num_out_eigval, aeroConverged)
!
! Use
!
  use globale_variablen, ONLY : shortoutFile
  use konstanten
  use fesimulation_typen
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
! Input
!
  type(fe_simulation)                                           :: fesim
  integer, intent(in)                                           :: num_out_eigval
  integer, intent(in)                                           :: sc
  double precision, dimension(num_out_eigval), intent(in)       :: eigenvalues
  double precision, dimension(fesim%internals%dim_dof,num_out_eigval), intent(in) :: eigenvectors
  logical, intent(in)                                           :: aeroConverged
!
! Internal
!
  integer                                                       :: ii
  integer, parameter                                            :: outfem = 55
  integer                                                       :: err_code=0
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
 
  if (textoutput .eq. .true.) then 
    write(*,*) 'Beulloesung fertig, LC', sc
    do ii = 1, num_out_eigval
      write(*,*) 'Subcase ', sc, ', Eigenwert ', ii, ': ', eigenvalues(ii)
    end do
  end if
  
  if (eigenvalues(1) .LT. 0.9d0) then
      fesim%internals%failed = .TRUE.
  end if
 
  if (fesim%ausgabe%outputVTK .eq. .true.) then
      call sol2_output_vtk(fesim, sc, eigenvalues, eigenvectors, num_out_eigval, aeroConverged)
  endif
  
  if (fesim%ausgabe%outputUser .eq. .true.) then
      call sol2_output_user(fesim, sc, eigenvalues, eigenvectors, num_out_eigval, aeroConverged)
  endif
  
  if (fesim%ausgabe%outputShort .eq. .true.) then
      OPEN(shortoutFile, file = 'output_FEM.txt', STATUS='OLD', POSITION='APPEND')
      write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', sc, ', Eigenwert             :  ', eigenvalues(1)
      CLOSE(shortoutFile)
  endif
  
  if ((fesim%ausgabe%outputHyMoWi .eq. .true.) .and. (aeroConverged .eq. .true.)) then
      OPEN(outfem, file = 'output_FEM_HyMoWi.txt', STATUS='OLD', POSITION='APPEND')
      write(outfem,'(A11,E25.18)') 'Eigenwert: ', eigenvalues(1)
      CLOSE(outfem)
  endif
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                       'sol2_output_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE sol2_output_control
