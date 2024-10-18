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
!> Create transformation matrix for strain vectors and orthotropic material
!> stiffness matrices out of a cosine matrix
!
!> @details
!> Creates the transformation matrix needed to transform a strain vector
!> or a orthotrop material stiffness matrix C to a coordinate system defined
!> by the cosine matrix E;
!> Strain vector or material stiffness matrix must correspond to the ordering
!> of the strain vector e as followed:
!> e = [e11,e22,e33,gamma23,gamma13,gamma12]^T
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_transformation_matrix.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_transformation_matrix(E,T)
! =================================================================================================
! use
!
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
double precision, intent(in)                  :: E(3,3) !< Cosine matrix for transformation
!
! Output
!
double precision, dimension(6,6), intent(out) :: T !< Transformation matrix
!
! Internal
!
integer                                       :: err_code=0
!
! =================================================================================================
! Initialize
 T = 0.d0
!
! =================================================================================================
!
! Calculation

 ! Set up full transformation matrix corresponding to the
 ! following ordering of the strain-vector e:
 ! e = [e11,e22,e33,gamma23,gamma13,gamma12]^T
 T(1,1) = E(1,1)**2
 T(1,2) = E(1,2)**2
 T(1,3) = E(1,3)**2
 T(1,4) = E(1,2)*E(1,3)
 T(1,5) = E(1,1)*E(1,3)
 T(1,6) = E(1,1)*E(1,2)
 T(2,1) = E(2,1)**2
 T(2,2) = E(2,2)**2
 T(2,3) = E(2,3)**2
 T(2,4) = E(2,2)*E(2,3)
 T(2,5) = E(2,1)*E(2,3)
 T(2,6) = E(2,1)*E(2,2)
 T(3,1) = E(3,1)**2
 T(3,2) = E(3,2)**2
 T(3,3) = E(3,3)**2
 T(3,4) = E(3,2)*E(3,3)
 T(3,5) = E(3,1)*E(3,3)
 T(3,6) = E(3,1)*E(3,2)
 T(4,1) = 2.d0*E(2,1)*E(3,1)
 T(4,2) = 2.d0*E(2,2)*E(3,2)
 T(4,3) = 2.d0*E(2,3)*E(3,3)
 T(4,4) = E(2,2)*E(3,3) + E(2,3)*E(3,2)
 T(4,5) = E(2,1)*E(3,3) + E(2,3)*E(3,1)
 T(4,6) = E(2,1)*E(3,2) + E(2,2)*E(3,1)
 T(5,1) = 2.d0*E(1,1)*E(3,1)
 T(5,2) = 2.d0*E(1,2)*E(3,2)
 T(5,3) = 2.d0*E(1,3)*E(3,3)
 T(5,4) = E(1,2)*E(3,3) + E(1,3)*E(3,2)
 T(5,5) = E(1,1)*E(3,3) + E(1,3)*E(3,1)
 T(5,6) = E(1,1)*E(3,2) + E(1,2)*E(3,1)
 T(6,1) = 2.d0*E(1,1)*E(2,1)
 T(6,2) = 2.d0*E(1,2)*E(2,2)
 T(6,3) = 2.d0*E(1,3)*E(2,3)
 T(6,4) = E(1,2)*E(2,3) + E(1,3)*E(2,2)
 T(6,5) = E(1,1)*E(2,3) + E(1,3)*E(2,1)
 T(6,6) = E(1,1)*E(2,2) + E(1,2)*E(2,1)
!
! =================================================================================================
!
! Error handling
!
9999 continue
!
if (err_code /= 0) then
!   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_transformation_matrix'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!   
end if
!
return
!
end subroutine lsolid20_transformation_matrix
