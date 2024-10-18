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
!> Calculates the nodal coordinate systems for a 8-node shell element quad8
!
!> @details
!> Calculation of the nodal coordinate systems for a 8-node shell element quad8
!
!> @author
!> Andreas Hauffe, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 23.05.2017
!
!> $Id: quad8_normals.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_normals(node_coords,normalVec,nodeTransMat,nodeTransMat2)
! =================================================================================================
! use
!  
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
!
! Output
!
  double precision, dimension(3,8), intent(out)     :: normalVec          !< Normalenvektoren an den acht Punkten
  double precision, dimension(8,3,3), intent(out)   :: nodeTransMat       !< Umrechnungsmatrix der Rotationen in Anteile der Verschiebung
  double precision, dimension(3,3,8), intent(out)   :: nodeTransMat2      !< lokales Koordinatensystem an jedem Knoten (Transformationmatrix lokal -> global)
!
! Internal
!
  double precision, dimension(8)                    :: dNidxi,dNideta
  double precision                                  :: dxdxi, dxdeta, dydxi, dydeta, dzdxi, dzdeta
  
  double precision, dimension(8)                    :: xi_n, eta_n
  
  double precision, dimension(3)                    :: x_loc, y_loc, z_loc      ! lokale Koordinatenachsen
  
  integer                                           :: ii
  
  integer                                           :: err_code=0
  
  ! Dieses Epsilon ist notwendig, um eine vollständige Übereinstimmung der Elementsteifigkeitsmatrix zwischen diesem Element und dem Shell 91 von Ansys zu erhalten
  double precision, parameter                       :: eps_natKoord = 0.001d0
  
  xi_n(1) = -1.d0+eps_natKoord; eta_n(1) = -1.d0+eps_natKoord
  xi_n(2) =  1.d0-eps_natKoord; eta_n(2) = -1.d0+eps_natKoord
  xi_n(3) =  1.d0-eps_natKoord; eta_n(3) =  1.d0-eps_natKoord
  xi_n(4) = -1.d0+eps_natKoord; eta_n(4) =  1.d0-eps_natKoord
  xi_n(5) =  0.d0             ; eta_n(5) = -1.d0+eps_natKoord
  xi_n(6) =  1.d0-eps_natKoord; eta_n(6) =  0.d0
  xi_n(7) =  0.d0             ; eta_n(7) =  1.d0-eps_natKoord
  xi_n(8) = -1.d0+eps_natKoord; eta_n(8) =  0.d0

  !*****************************************************
  ! Berechnung des Normalenvektors in den Knotenpunkten
  !*****************************************************
  do ii = 1,8
  
    call quad8_ansatzfunction_xieta(xi_n(ii), eta_n(ii), dNidxi=dNidxi, dNideta=dNideta)
    
    dxdxi  = SUM(dNidxi(1:8) *node_coords(1:8,1))
    dxdeta = SUM(dNideta(1:8)*node_coords(1:8,1))
    
    dydxi  = SUM(dNidxi(1:8) *node_coords(1:8,2))
    dydeta = SUM(dNideta(1:8)*node_coords(1:8,2))
    
    dzdxi  = SUM(dNidxi(1:8) *node_coords(1:8,3))
    dzdeta = SUM(dNideta(1:8)*node_coords(1:8,3))
     
    ! Berechnung der Normalen auf die Elementmittelsfläche
    normalVec(1,ii) = dydxi*dzdeta - dzdxi*dydeta
    normalVec(2,ii) = dzdxi*dxdeta - dxdxi*dzdeta
    normalVec(3,ii) = dxdxi*dydeta - dydxi*dxdeta
    
    normalVec(:,ii) = normalVec(:,ii) / SQRT(DOT_PRODUCT(normalVec(:,ii),normalVec(:,ii)))
    
    ! lokales Koordinatensystem berechnen
    
    x_loc(1) = dxdxi
    x_loc(2) = dydxi
    x_loc(3) = dzdxi
    
    x_loc = x_loc / SQRT(DOT_PRODUCT(x_loc,x_loc))
    
    ! Berechnung der Normalen auf die Elementmittelsfläche am Integrationspunkt
    z_loc(1:3) = normalVec(1:3,ii)
    
    ! cross_prod(z_vec, x_vec)
    y_loc(1) = z_loc(2)*x_loc(3) - x_loc(2)*z_loc(3)
    y_loc(2) = z_loc(3)*x_loc(1) - x_loc(3)*z_loc(1)
    y_loc(3) = z_loc(1)*x_loc(2) - x_loc(1)*z_loc(2)
    y_loc = y_loc / SQRT(DOT_PRODUCT(y_loc,y_loc))
    
    ! Umrechnung der lokalen Rotationen in Anteile der Verschiebung
    nodeTransMat(ii,1:3,1) = -y_loc(1:3)
    nodeTransMat(ii,1:3,2) =  x_loc(1:3)
    nodeTransMat(ii,1:3,3) = 0.d0
    
    ! Transformationsmatrix, um die Rotationen vom lokalen System in das globale zu transformieren
    nodeTransMat2(1:3,1,ii) = x_loc(1:3)
    nodeTransMat2(1:3,2,ii) = y_loc(1:3)
    nodeTransMat2(1:3,3,ii) = z_loc(1:3)
    
    nodeTransMat(ii,1:3,1:3) = matmul(nodeTransMat(ii,1:3,1:3),transpose(nodeTransMat2(1:3,1:3,ii)))
    
  end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_normals'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_normals
