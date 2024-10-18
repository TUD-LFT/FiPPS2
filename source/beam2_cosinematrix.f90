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
!> Subroutine for transformation matrix for beam2-element
!
!> @details 
!> Subroutines computes the transformation matrix from local to global
!> coordinates (or the other way around) for a 2-node beam element
!> Additionally he element length is computed
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 22.06.2010
!
!> $Id: beam2_rotation.f90 362 2018-08-07 09:18:14Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 362 $
!> $Date: 2018-08-07 11:18:14 +0200 (Di, 07. Aug 2018) $
! 
! =================================================================================================
subroutine beam2_cosinematrix (beam2_node_coords,vv,BC,le)

use konstanten
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
double precision, dimension (2,3), intent(in)    :: beam2_node_coords   !< nodal coordinates of element nodes in global coordinates
double precision, dimension(3), intent(in)       :: vv                  !< vector in x-y-plane of local element coordinate system
!
! Output
!
double precision, dimension (3,3), intent(out)   :: BC                  !< cosine matrix
double precision, intent(out)                    :: le                  !< element length
!
! Input + Output
!

!
! Internal
!
double precision, dimension(3)                   :: AB,e1,e2,e3,vec3
integer                                          :: variante_TRMat
double precision                                 :: dd,S1,S2,S3,C1,C2,C3,lxy

integer                                          :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
! =================================================================================================
!
! Calculation
!
  ! calculation of vector along x-axis of beam
  
  AB = beam2_node_coords(2,1:3)-beam2_node_coords(1,1:3)
  
  ! element length

  le=sqrt(dot_product(AB,AB))
  
  ! transformation matrix
  
  variante_TRMat=1
  if (SUM(vv) .EQ. 0.d0) variante_TRMat=2
  
  select case (variante_TRMat)
    
    case(1) ! NASTRAN
      
      ! unit-vector in x-direction
      
      e1 = AB/le
      
      ! unit vector in z-direction

      vec3(1) = AB(2)*vv(3) - AB(3)*vv(2)
      vec3(2) = AB(3)*vv(1) - AB(1)*vv(3)
      vec3(3) = AB(1)*vv(2) - AB(2)*vv(1)

      e3=vec3/sqrt(dot_product(vec3,vec3))
      
      ! unit vector in y-direction

      e2(1) = -(e1(2)*e3(3) - e1(3)*e3(2))
      e2(2) = -(e1(3)*e3(1) - e1(1)*e3(3))
      e2(3) = -(e1(1)*e3(2) - e1(2)*e3(1))

      !
      ! transformation submatrix
      !

      BC(1:3,1) = e1
      BC(1:3,2) = e2
      BC(1:3,3) = e3
      
    case(2) ! ANSYS
    
      ! advantage of this method is, that it does not need an auxiliary vector vv
      
      dd = 0.0001d0*le
      
      lxy = sqrt((beam2_node_coords(2,1)-beam2_node_coords(1,1))**2.D0+(beam2_node_coords(2,2)-beam2_node_coords(1,2))**2.D0)
      
      if (lxy .GT. dd) then
        
       S1 = (beam2_node_coords(2,2)-beam2_node_coords(1,2))/lxy
       C1 = (beam2_node_coords(2,1)-beam2_node_coords(1,1))/lxy

      else
        
       S1 = 0.D0
       C1 = 1.D0
        
      end if
      
      S2 = (beam2_node_coords(2,3)-beam2_node_coords(1,3))/le
      C2 = lxy/le
      
      S3 = 0.D0
      C3 = 1.D0
      
      BC(1,1) =  C1*C2
      BC(1,2) =  S1*C2
      BC(1,3) =  S2
      BC(2,1) = -C1*S2*S3-S1*C3
      BC(2,2) = -S1*S2*S3+C1*C3
      BC(2,3) =  S3*C2
      BC(3,1) = -C1*S2*C3-S1*S3
      BC(3,2) = -S1*S2*C3-C1*S3
      BC(3,3) =  C3*C2
      
      ! Trafo-Matrix given in ANSYS-Help sec. 14.4.4 (Element Library -> BEAM4 3D Elastic Beam -> Local to Global Conversion
      ! is for transformation from global to local
      ! for transformation from local to global the transposed matrix has to be used
      
!      if (mode == 'lg') then
        
      BC = transpose(BC)

!      end if

  end select

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_cosinematrix'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine beam2_cosinematrix
