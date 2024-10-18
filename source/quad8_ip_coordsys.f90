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
!> Calculates the transformation matrix from the element coordinate system at the integration point
!> to the global system
!
!> @details
!> Calculates the transformation matrix from the element coordinate system at the integration point
!> to the global system with respect to the element x direction
!
!> @author
!> Andreas Hauffe, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 23.05.2017
!
!> $Id: quad8_ip_coordsys.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_ip_coordsys(node_coords, xi, eta, xElem, transMat)
! =================================================================================================
! use
!
  use integration_schemes
  
  implicit none
!
! =================================================================================================
!
! Interface
!
  interface 
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
      double precision, intent(in)                              :: xi, eta
      
      double precision, dimension(8), optional, intent(out)     :: Ni
      double precision, dimension(8), optional, intent(out)     :: dNidxi, dNideta
      double precision, dimension(8), optional, intent(out)     :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
  end interface
!
! Input
!
  double precision, intent(in)                      :: node_coords(8,3)   !< (8,3)-Array containing global coordinates of element nodes
  double precision                                  :: xi                 !< xi coordinate of the integration point
  double precision                                  :: eta                !< eta coordinate of the integration point
  double precision, dimension(3)                    :: xElem              !< x direction of the element coordinate system
!
! Output
!
  double precision, dimension(3,3), intent(out)     :: transMat           !< 3x3 transformation matrix from local integration point element coordinate system to global system
!
! Internal
!
  double precision, dimension(8)                    :: dNidxi,dNideta
  double precision                                  :: dxdxi, dxdeta, dydxi, dydeta, dzdxi, dzdeta
  
  double precision, dimension(3)                    :: x_loc, y_loc, z_loc      ! lokale Koordinatenachsen
  
  integer                                           :: err_code=0
  
  call quad8_ansatzfunction_xieta(xi, eta, dNidxi=dNidxi, dNideta=dNideta)
  
  dxdxi  = SUM(dNidxi(1:8) *node_coords(1:8,1))
  dxdeta = SUM(dNideta(1:8)*node_coords(1:8,1))
  
  dydxi  = SUM(dNidxi(1:8) *node_coords(1:8,2))
  dydeta = SUM(dNideta(1:8)*node_coords(1:8,2))
  
  dzdxi  = SUM(dNidxi(1:8) *node_coords(1:8,3))
  dzdeta = SUM(dNideta(1:8)*node_coords(1:8,3))
   
  ! Berechnung der Normalen auf die Elementmittelsfläche
  z_loc(1) = dydxi*dzdeta - dzdxi*dydeta
  z_loc(2) = dzdxi*dxdeta - dxdxi*dzdeta
  z_loc(3) = dxdxi*dydeta - dydxi*dxdeta
  z_loc(:) = z_loc(:) / SQRT(DOT_PRODUCT(z_loc(:),z_loc(:)))
  
  ! cross_prod(z_vec, x_vec)
  y_loc(1) = z_loc(2)*xElem(3) - xElem(2)*z_loc(3)
  y_loc(2) = z_loc(3)*xElem(1) - xElem(3)*z_loc(1)
  y_loc(3) = z_loc(1)*xElem(2) - xElem(1)*z_loc(2)
  y_loc = y_loc / SQRT(DOT_PRODUCT(y_loc,y_loc))
  
  ! cross_prod(y_vec, z_vec)
  x_loc(1) = y_loc(2)*z_loc(3) - z_loc(2)*y_loc(3)
  x_loc(2) = y_loc(3)*z_loc(1) - z_loc(3)*y_loc(1)
  x_loc(3) = y_loc(1)*z_loc(2) - z_loc(1)*y_loc(2)
  x_loc = x_loc / SQRT(DOT_PRODUCT(x_loc,x_loc))

  ! Tranformationsmatrix - Element- bzw. Materialkoordinatensystem in globales System
  transMat(1:3,1) = x_loc(1:3)
  transMat(1:3,2) = y_loc(1:3)
  transMat(1:3,3) = z_loc(1:3)

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_ip_coordsys'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_ip_coordsys
