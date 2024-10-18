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
!> Set up cosine matrix for element coordinate system using a given jacobi-matrix
!
!> @details
!> Cosine matrix for element coordinate system is calculated with a given jacobi
!> matrix at a specific point inside the element; if a local coordinate system is
!> given (cid > 0), the x-axis of the element coordinate system is aligned along
!> the x-axis of the local coordinate system; alternatively, cid = -1 makes use of
!> the x-axis of the global coordinate system; if cid = 0, the orientation of the
!> element coordinate system is only specified by the element's orientation, so that
!> the x-axis of the element coordinate system is aligned along the xi-coordinate;
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> @author Andreas Hauffe, TU Dresden, Zeilen und Spalten von E(:,:) vertauscht, um das Erzeugen von temporären Arrays für die Übergabe in cross_product3 zu verhindern
!
!> $Id: cosinematrix_esys.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine cosinematrix_esys(fesim, J, cid, E)
! =================================================================================================
! use
 use konstanten
 use mat_func
 use fesimulation_typen
!
!
! =================================================================================================
!
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
! Input
!
type(fe_simulation), intent(in)               :: fesim !< FE-Simulationsspeicher
integer, intent(in)                           :: cid !< ID of local coordinate system used for orientation of the x-axis of the element coordinate system
double precision, intent(in), dimension(3,3)  :: J !< Jacobi matrix
!
! Output
!
double precision, intent(out), dimension(3,3) :: E !< Cosine matrix of the element coordinate system
!
! Internal
!
integer                          :: err_code=0
double precision                 :: n1, n2
double precision, dimension(3)   :: e_orient
double precision, dimension(3,3) :: T
!
! =================================================================================================
!
! Calculation
!
 ! Set up cosine matrix for transforming material stiffness matrix from element to
 ! global coordinate system;
 ! Create local orthonormal coordinate system (element coordinate system)
 ! with first unit vector aligned along xi-direction, third unit vector points
 ! in thickness direction
 ! First unit vector T(:,1)
 T(:,1) = J(1,:)
 T(:,1) = T(:,1) / sqrt(dot_product(T(:,1),T(:,1)))
 ! Temporary vector
 T(:,2) = J(2,:)
 ! Third unit vector T(:,3)
 T(:,3) = cross_product3(T(:,1),T(:,2))
 T(:,3) = T(:,3) / sqrt(dot_product(T(:,3),T(:,3)))
 ! Second unit vector T(:,2)
 T(:,2) = cross_product3(T(:,3),T(:,1))
 T(:,2) = T(:,2) / sqrt(dot_product(T(:,2),T(:,2)))
 ! Check, if user defined orientation shall be used
 if (cid .NE. 0) then
   if (cid .EQ. -1) then
     ! Get x-axis-vector of global coosy
     e_orient(:) = (/ 1.d0, 0.d0, 0.d0 /)
   else
     ! Check if input is valid
     if ((cid .LT. 1) .OR. (cid .GT. size(fesim%koordinatensysteme%coords,1))) then
       write(*,*) 'Error orientating element coordinate system'
       write(*,*) 'Invalid input for orientation coordinate system'
       write(*,*) 'ID of orientation coordinate system must specifiy'
       write(*,*) 'a defined coordinate system'
       write(*,*) 'Alternatively, use cid = -1 for orientation along'
       write(*,*) 'global coordinate system or cid = 0 for standard'
       write(*,*) 'orientation'
       err_code = 2
       goto 9999
     else
       ! Get x-axis-vector of local coosy
       e_orient(:) = fesim%koordinatensysteme%coords(cid)%transMat(1:3,1)
     end if
   end if
   ! Check if angle between x-axis-vector and element's thickness direction is
   ! greater than 5 degrees, otherwise throw error
   n1 = sqrt(dot_product(e_orient(:),e_orient(:)))
   n2 = sqrt(dot_product(T(:,3),T(:,3)))
   if (abs(dot_product(e_orient,T(:,3))/(n1*n2)) .GE. cos(5.d0*pi/180.d0)) then
     write(*,*) 'Error computing local orthonormal element coordinate system'
     write(*,*) "Angle between vector prescribed as direction for element's"
     write(*,*) 'x-axis and vector along thickness direction is less than 5'
     write(*,*) 'degree. Please check the orientation of the vector.'
     err_code = 1
     goto 9999
   end if
   ! new first unit vector as projection of e_orient in plane of T(:,1) and T(:,2)
   T(:,1) = dot_product(e_orient(:),T(:,1))*T(:,1) + dot_product(e_orient(:),T(:,2))*T(:,2)
   T(:,1) = T(:,1) / sqrt(dot_product(T(:,1),T(:,1)))
   ! new second unit vector
   T(:,2) = cross_product3(T(:,3),T(:,1))
   T(:,2) = T(:,2) / sqrt(dot_product(T(:,2),T(:,2)))
   ! third unit vector remains unchanged
 end if
 
 E(1,:) = T(:,1)
 E(2,:) = T(:,2)
 E(3,:) = T(:,3)
 
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'cosinematrix_esys'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine cosinematrix_esys
