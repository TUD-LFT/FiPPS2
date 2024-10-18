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
subroutine quad8_iso_bmat(xi, eta, zeta, node_coords, thickness, normalVec, nodeTransMat, bmat, detJac, transMat)

  implicit none
!
! =================================================================================================
!
! Interface
!
  interface 
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
      double precision, intent(in)				    :: xi, eta
      
      double precision, dimension(8), optional, intent(out)	    :: Ni
      double precision, dimension(8), optional, intent(out)	    :: dNidxi, dNideta
      double precision, dimension(8), optional, intent(out)	    :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
  end interface
  
    
  double precision, intent(in)                    :: xi, eta, zeta
  double precision, dimension(8,3), intent(in) 	  :: node_coords
  double precision, intent(in)			  :: thickness
  double precision, dimension(8,3,3), intent(in)  :: nodeTransMat		! lokales Koordinatensystem an jedem Knoten
  
  double precision, dimension(48,6), intent(out)  :: bMat
  double precision, intent(out)                   :: detJac

  double precision, dimension(8)		  :: Ni, dNidxi,dNideta
  double precision                                :: dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta, dzdxi, dzdeta, dzdzeta
  double precision, dimension(3,8)      	  :: normalVec 			! Normalenvektoren an den acht Punkten
  double precision, dimension(3,3)                :: jac, jacInv		! Jakobi-Matrix und inverse Jakobi-Matrix an den Integrationspunkten
  double precision, dimension(3)                  :: x_loc, y_loc, z_loc	! lokale Koordinatenachsen
  double precision, dimension(3,3), intent(out)   :: transMat
  double precision, dimension(8,3)                :: n1, n2

  logical	                                  :: ok_flag

  call quad8_ansatzfunction_xieta(xi, eta, Ni=Ni, dNidxi=dNidxi, dNideta=dNideta)
  
  ! Einträge der Jakobimatrix bestimmen
  
  dxdxi   = 0.D0
  dxdeta  = 0.D0
  dxdzeta = 0.D0
  dydxi   = 0.D0
  dydeta  = 0.D0
  dydzeta = 0.D0
  dzdxi   = 0.D0
  dzdeta  = 0.D0
  dzdzeta = 0.D0
  
  dxdxi   = SUM(dNidxi(1:8) * (node_coords(1:8,1) + thickness/2.d0*zeta*normalVec(1,1:8)))
  dxdeta  = SUM(dNideta(1:8)* (node_coords(1:8,1) + thickness/2.d0*zeta*normalVec(1,1:8)))
  dxdzeta = SUM(Ni(1:8)*thickness/2.d0*normalVec(1,1:8))
  
  dydxi   = SUM(dNidxi(1:8) * (node_coords(1:8,2) + thickness/2.d0*zeta*normalVec(2,1:8)))
  dydeta  = SUM(dNideta(1:8)* (node_coords(1:8,2) + thickness/2.d0*zeta*normalVec(2,1:8)))
  dydzeta = SUM(Ni(1:8)*thickness/2.d0*normalVec(2,1:8))
  
  dzdxi   = SUM(dNidxi(1:8) * (node_coords(1:8,3) + thickness/2.d0*zeta*normalVec(3,1:8)))
  dzdeta  = SUM(dNideta(1:8)* (node_coords(1:8,3) + thickness/2.d0*zeta*normalVec(3,1:8)))
  dzdzeta = SUM(Ni(1:8)*thickness/2.d0*normalVec(3,1:8))
  
  jac(1,1) =   dxdxi; jac(1,2) =   dydxi; jac(1,3) =   dzdxi;
  jac(2,1) =  dxdeta; jac(2,2) =  dydeta; jac(2,3) =  dzdeta;
  jac(3,1) = dxdzeta; jac(3,2) = dydzeta; jac(3,3) = dzdzeta;
    
  call m33inv(jac, jacInv, detJac, ok_flag)
  
  ! lokales Koordinatensystem berechnen
  
  x_loc(1) = dxdxi
  x_loc(2) = dydxi
  x_loc(3) = dzdxi
  
  x_loc = x_loc / SQRT(DOT_PRODUCT(x_loc,x_loc))
  
  ! Berechnung der Normalen auf die Elementmittelsfläche am Integrationspunkt
  z_loc(1) = dydxi*dzdeta - dzdxi*dydeta
  z_loc(2) = dzdxi*dxdeta - dxdxi*dzdeta
  z_loc(3) = dxdxi*dydeta - dydxi*dxdeta
  
  z_loc = z_loc / SQRT(DOT_PRODUCT(z_loc,z_loc))
  
  ! cross_prod(z_vec, x_vec) (sollte normiert sein, da x und z normiert wurden)
  y_loc(1) = z_loc(2)*x_loc(3) - x_loc(2)*z_loc(3)
  y_loc(2) = z_loc(3)*x_loc(1) - x_loc(3)*z_loc(1)
  y_loc(3) = z_loc(1)*x_loc(2) - x_loc(1)*z_loc(2)

  transMat(1:3,1) = x_loc(1:3)
  transMat(1:3,2) = y_loc(1:3)
  transMat(1:3,3) = z_loc(1:3)
  
  jacInv = matmul(transpose(transMat), jacInv)
  
  n1(1:8,1) = dNidxi(1:8)
  n1(1:8,2) = dNideta(1:8)
  n1(1:8,3) = 0.d0
  
  n2(1:8,1) = 0.5d0 * thickness * dNidxi(1:8)  * zeta
  n2(1:8,2) = 0.5d0 * thickness * dNideta(1:8) * zeta
  n2(1:8,3) = 0.5d0 * thickness * Ni(1:8)
  
  ! Durch MatLab (quad8_iso_emat.m) mehr oder weniger automatisch generierter Code, um Operationen mit Null zu vermeiden
  
  bMat( 1: 8, 1) = transMat(1, 1)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
  bMat( 9:16, 1) = transMat(2, 1)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
  bMat(17:24, 1) = transMat(3, 1)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
  bMat(25:32, 1) = transMat(1, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
  bMat(33:40, 1) = transMat(1, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
  bMat(41:48, 1) = transMat(1, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
  bMat( 1: 8, 2) = transMat(1, 2)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
  bMat( 9:16, 2) = transMat(2, 2)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
  bMat(17:24, 2) = transMat(3, 2)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
  bMat(25:32, 2) = transMat(1, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
  bMat(33:40, 2) = transMat(1, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
  bMat(41:48, 2) = transMat(1, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
  bMat( 1: 8, 3) = transMat(1, 3)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat( 9:16, 3) = transMat(2, 3)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat(17:24, 3) = transMat(3, 3)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat(25:32, 3) = transMat(1, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
  bMat(33:40, 3) = transMat(1, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
  bMat(41:48, 3) = transMat(1, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
  bMat( 1: 8, 4) = transMat(1, 2)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2)) + transMat(1, 1)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
  bMat( 9:16, 4) = transMat(2, 2)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2)) + transMat(2, 1)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
  bMat(17:24, 4) = transMat(3, 2)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2)) + transMat(3, 1)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
  bMat(25:32, 4) = transMat(1, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(1, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(2, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1) + transMat(3, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
  bMat(33:40, 4) = transMat(1, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(1, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(2, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2) + transMat(3, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
  bMat(41:48, 4) = transMat(1, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(1, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(2, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3) + transMat(3, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
  bMat( 1: 8, 5) = transMat(1, 3)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2)) + transMat(1, 2)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat( 9:16, 5) = transMat(2, 3)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2)) + transMat(2, 2)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat(17:24, 5) = transMat(3, 3)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2)) + transMat(3, 2)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat(25:32, 5) = transMat(1, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(1, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(2, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1) + transMat(3, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
  bMat(33:40, 5) = transMat(1, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(1, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(2, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2) + transMat(3, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
  bMat(41:48, 5) = transMat(1, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(1, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(2, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3) + transMat(3, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
  bMat( 1: 8, 6) = transMat(1, 3)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2)) + transMat(1, 1)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat( 9:16, 6) = transMat(2, 3)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2)) + transMat(2, 1)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat(17:24, 6) = transMat(3, 3)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2)) + transMat(3, 1)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
  bMat(25:32, 6) = transMat(1, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(1, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(3, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1) + transMat(2, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
  bMat(33:40, 6) = transMat(1, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(1, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(3, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2) + transMat(2, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
  bMat(41:48, 6) = transMat(1, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(1, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(3, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3) + transMat(2, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    
  contains
!***********************************************************************************************************************************
!  Quelle: http://www.davidgsimpson.com/software/m33inv_f90.txt
!
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************

    SUBROUTINE M33INV (A, AINV, DET, OK_FLAG)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      DOUBLE PRECISION, INTENT(OUT) :: DET
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

    END SUBROUTINE M33INV
    
end subroutine quad8_iso_bmat
