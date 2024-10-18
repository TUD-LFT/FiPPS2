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
!> Calculates the area of the reference surface of a quad8 finite
!> shell element
!
!> @details
!               
!> @author   Andreas Hauffe, TU Dresden, wiss. Mitarbeiter
!
!> $Id: quad8_area.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_area(node_coords,area)

  use integration_schemes
  use konstanten

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
  
  double precision, dimension(8,3), intent(in)    :: node_coords  !< (8,3)-Array containing the coordinate of each node
  double precision, intent(out)                   :: area         !< Area of the reference surface

  double precision, dimension(8)                  :: dNidxi,dNideta
  double precision                                :: dxdxi, dxdeta, dydxi, dydeta, dzdxi, dzdeta
  
  double precision, dimension(:), allocatable     :: xi_vec, eta_vec, w_xi, w_eta
  
  double precision                                :: Exieta, Fxieta, Gxieta
  
  integer                                         :: ii
  integer                                         :: numgausspoints
  

  !*****************************************************
  ! Integration
  !*****************************************************
  
  numgausspoints = nip_quad8*nip_quad8
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))
  
  call integration_points_2d(nip_quad8, "Gauss",  xi_vec, w_xi, &
                             nip_quad8, "Gauss", eta_vec, w_eta)
  
  ! Berechnung der Fläche des Elements
  ! Sollte eventuell mit einer 3x3 Integration durchgeführt werden
  
  area = 0.d0
  
  do ii = 1,nip_quad8*nip_quad8
  
    call quad8_ansatzfunction_xieta(xi_vec(ii), eta_vec(ii), dNidxi=dNidxi, dNideta=dNideta)
    
    dxdxi  = 0.D0
    dxdeta = 0.D0
    dydxi  = 0.D0
    dydeta = 0.D0
    dzdxi  = 0.D0
    dzdeta = 0.D0
    
    dxdxi  = SUM(dNidxi(1:8) *node_coords(1:8,1))
    dxdeta = SUM(dNideta(1:8)*node_coords(1:8,1))
    
    dydxi  = SUM(dNidxi(1:8) *node_coords(1:8,2))
    dydeta = SUM(dNideta(1:8)*node_coords(1:8,2))
    
    dzdxi  = SUM(dNidxi(1:8) *node_coords(1:8,3))
    dzdeta = SUM(dNideta(1:8)*node_coords(1:8,3))
      
    Exieta = dxdxi**2 + dydxi**2 + dzdxi**2
    Fxieta = dxdxi*dxdeta + dydxi*dydeta + dzdxi*dzdeta
    Gxieta = dxdeta**2 + dydeta**2 + dzdeta**2
    
    area = area + w_xi(ii)*w_eta(ii)*SQRT(Exieta*Gxieta-Fxieta**2)
  
  end do

  deallocate(xi_vec, eta_vec)
  deallocate(  w_xi,   w_eta)

end subroutine quad8_area