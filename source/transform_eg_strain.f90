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
!> Transform strain vector from element to global coordinate system or reverse
!> using a given jacobi-matrix
!
!> @details
!> Strain vector is transformed from element to global coordinate system or reverse
!> using the jacobi matrix which must be given at the point of the element, where
!> transformation shall take place; ordering of the strain vector must be as followed:
!> e = [e11, e22, e33, e23, e13, e12]^T
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: transform_eg_strain.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine transform_eg_strain(fesim, eps_orig, J, cid, eps_targ, dir)
! =================================================================================================
! use
!
 use fesimulation_typen
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
type(fe_simulation), intent(in)               :: fesim
double precision, intent(in), dimension(6)    :: eps_orig !< Strain vector with respect to the original coordinate system
double precision, intent(in), dimension(3,3)  :: J !< Jacobi matrix at the point where transformation shall take place
integer, intent(in)                           :: cid !< ID of coordinate system for orientation of element coosy
character(2), intent(in)                      :: dir !< Specifies the direction of the transformation: dir = 'eg': element -> global  coordinate system, dir = 'ge': global  -> element coordinate system
!
! Output
!
double precision, intent(out), dimension(6)   :: eps_targ !< Strain vector with respect to the target coordinate system
!
! Internal
!
integer                          :: err_code=0
double precision, dimension(3,3) :: E
double precision, dimension(6,6) :: T
!
! =================================================================================================
!
! Calculation
!
 ! Get cosine matrix for transforming material stiffness matrix from element to
 ! global coordinate system;
 call cosinematrix_esys(fesim,J,cid,E)

 if (dir .EQ. 'eg') then
   ! Transformation of strain needs transposed E-matrix
   E = transpose(E)
 else if (dir .EQ. 'ge') then
   ! do nothing
 else
   write(*,*) 'wrong specifier for direction of coordinate'
   write(*,*) 'transformation is chosen'
   write(*,*) 'valid specifiers are:'
   write(*,*) 'dir = ''eg'': element -> global  coordinate system'
   write(*,*) 'dir = ''ge'': global  -> element coordinate system'
   err_code = 2
   goto 9999
 end if
 
 ! Get Transformation matrix T
 call lsolid20_transformation_matrix(E,T)
 
 ! Transform strain
 eps_targ = matmul(T,eps_orig)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'transform_eg_strain'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine transform_eg_strain
