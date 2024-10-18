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
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 22.06.2010
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter 24.06.2010
!>                                  -> Auswahl Format ANSYS oder FiPPS hinzugefügt
!
!> $Id: write_quad_matrix.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine write_quad_matrix(mat, size, name, way)

  implicit none
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
  integer                                               :: ii, jj, kk
  integer, intent(in)                                   :: size
  
  double precision, dimension(size,size), intent(in)    :: mat
  
  character(20), intent(in)                             :: name
  character(5), intent(in)                              :: way  !< Darstellungsart - kann die Werte 'FiPPS','ANSYS' oder 'beide' annehmen
  character(20)                                         :: form
  
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
  if (way == 'FiPPS') then
  
    write(*,*) name
    WRITE(form,'(A1,I10.10,A4)') '(', size,'I14)'
    write(*,form) (jj, jj = 1,size)
    WRITE(form,'(A4,I10.10,A6)') '(I8,', size,'E14.7)'
    do ii = 1,size
      write(*,form) ii,(mat(ii,jj), jj = 1,size)
    end do
  
  else if (way == 'ANSYS') then
    
    write(*,*) name
    !! Matrix, wie Ansys herausschreiben
    do ii = 1,size
      do jj = 1,size/6
        write(*,'(I4,2X,6(E14.7,X))') ii,(mat(ii,kk), kk = (jj-1)*6+1,jj*6)
      end do
    end do
    
  else if (way == 'beide') then
    
    write(*,*) name
    
    WRITE(form,'(A1,I10.10,A4)') '(', size,'I14)'
    write(*,form) (jj, jj = 1,size)
    WRITE(form,'(A4,I10.10,A6)') '(I8,', size,'E14.7)'
    do ii = 1,size
      write(*,form) ii,(mat(ii,jj), jj = 1,size)
    end do
    
    !! Matrix, wie Ansys herausschreiben
    do ii = 1,size
      do jj = 1,size/6
        write(*,'(I4,2X,6(E14.7,X))') ii,(mat(ii,kk), kk = (jj-1)*6+1,jj*6)
      end do
    end do
    
  else
    
    write(*,*) 'wrong input on parameter way'
    write(*,*) 'this parameter can only take the following values:'
    write(*,*) 'FiPPS, ANSYS, beide'
    err_code=1
    goto 9999
    
  end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'write_quad_matrix'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine write_quad_matrix