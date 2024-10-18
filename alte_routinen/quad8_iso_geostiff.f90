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
subroutine quad8_iso_geostiff(disp_q8,E,nue,thickness,area,node_coords,quad8_Kg)

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
  
  double precision, intent(in)			  :: E,nue,thickness, area
  double precision, dimension(48), intent(in)     :: disp_q8
  double precision, dimension(8,3), intent(in) 	  :: node_coords
  double precision, dimension(48,48), intent(out) :: quad8_Kg

  double precision, dimension(8)		  :: Ni, dNidxi,dNideta
  double precision                                :: dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta, dzdxi, dzdeta, dzdzeta
  double precision, dimension(3,8)      	  :: normalVec 			! Normalenvektoren an den acht Punkten
  double precision, dimension(3)                  :: x_loc, y_loc, z_loc	! lokale Koordinatenachsen
  double precision, dimension(8,3,3)              :: nodeTransMat		! lokales Koordinatensystem an jedem Knoten
  
  double precision, dimension(8)        	  :: xi_n, eta_n
  
  double precision, dimension(:), allocatable     :: xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta
  
  double precision, dimension(48,48)              :: quad8_Kgt			! temporäre Elementsteifigkeitsmatrix
  
  double precision, dimension(48)                 :: reorderedDispl
  
  double precision, dimension(6)                  :: Nc
  
  double precision, dimension(3,3)                :: lambda
  
  double precision, dimension(9,9)                :: Sm
  
  double precision, dimension(3,3)                :: jac, jacInv		! Jakobi-Matrix und inverse Jakobi-Matrix an den Integrationspunkten
  
  double precision, dimension(3,3)                :: transMat
  
  double precision, dimension(3,3,8)              :: nodeTransMat2
  
  double precision, dimension(8,3)                :: n1, n2
  
  double precision, dimension(48,9)               :: bt
  
  double precision, dimension(6,9)                :: Hmat
  
  double precision, dimension(6,48)               :: bMat
  
  double precision                                :: fac, k, detJac
  
  double precision, dimension(6,6)                :: Dmat
  
  integer, parameter                              :: inplaneintppd  = 2		! Integrationspunkte in der Elementebene
  integer, parameter                              :: outplaneintppd = 2		! Integrationspunkte in Dickenrichtung
  
  integer                                         :: ii, jj
  integer                                         :: numgausspoints
  
  integer, dimension(48), parameter               :: index = (/ 1, 7,13,19,25,31,37,43, 2, 8,14,20,26,32,38,44, 3, 9,15,21,27,33,39,45, 4,10,16,22,28,34,40,46, 5,11,17,23,29,35,41,47, 6,12,18,24,30,36,42,48/)
  
  logical	                                  :: ok_flag

  xi_n(1) = -1.d0; eta_n(1) = -1.d0
  xi_n(2) =  1.d0; eta_n(2) = -1.d0
  xi_n(3) =  1.d0; eta_n(3) =  1.d0
  xi_n(4) = -1.d0; eta_n(4) =  1.d0
  xi_n(5) =  0.d0; eta_n(5) = -1.d0
  xi_n(6) =  1.d0; eta_n(6) =  0.d0
  xi_n(7) =  0.d0; eta_n(7) =  1.d0
  xi_n(8) = -1.d0; eta_n(8) =  0.d0
  
  do jj = 1,6
    do ii = 1,8
      reorderedDispl((jj-1)*8+ii) = disp_q8((ii-1)*6+jj)
    end do
  end do
  
  !*****************************************************
  ! Berechnung des Normalenvektors in den Knotenpunkten
  !*****************************************************
  do ii = 1,8
  
    call quad8_ansatzfunction_xieta(xi_n(ii), eta_n(ii), dNidxi=dNidxi, dNideta=dNideta)
    
    dxdxi  = 0.D0
    dxdeta = 0.D0
    dydxi  = 0.D0
    dydeta = 0.D0
    dzdxi  = 0.D0
    dzdeta = 0.D0
    
    do jj=1,8
    
      dxdxi  = dxdxi +dNidxi(jj) *node_coords(jj,1)
      dxdeta = dxdeta+dNideta(jj)*node_coords(jj,1)
    
      dydxi  = dydxi +dNidxi(jj) *node_coords(jj,2)
      dydeta = dydeta+dNideta(jj)*node_coords(jj,2)
    
      dzdxi  = dzdxi +dNidxi(jj) *node_coords(jj,3)
      dzdeta = dzdeta+dNideta(jj)*node_coords(jj,3)
    
    end do
    
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
    
    ! cross_prod(z_vec, x_vec) (sollte normiert sein, da x und z normiert wurden)
    y_loc(1) = z_loc(2)*x_loc(3) - x_loc(2)*z_loc(3)
    y_loc(2) = z_loc(3)*x_loc(1) - x_loc(3)*z_loc(1)
    y_loc(3) = z_loc(1)*x_loc(2) - x_loc(1)*z_loc(2)

    nodeTransMat(ii,1:3,1) = -y_loc(1:3)
    nodeTransMat(ii,1:3,2) =  x_loc(1:3)
    nodeTransMat(ii,1:3,3) = 0.d0 
    
    nodeTransMat2(1:3,1,ii) = x_loc(1:3)
    nodeTransMat2(1:3,2,ii) = y_loc(1:3)
    nodeTransMat2(1:3,3,ii) = z_loc(1:3)
    
    nodeTransMat(ii,1:3,1:3) = matmul(nodeTransMat(ii,1:3,1:3),transpose(nodeTransMat2(1:3,1:3,ii)))
    
  end do

  !*****************************************************
  ! Integration
  !*****************************************************
  
  numgausspoints = inplaneintppd*inplaneintppd*outplaneintppd                                           ! total number of Gauss Points
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints), zeta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints),   w_zeta(numgausspoints))
  
  call gauss_integration_values_3d (inplaneintppd, outplaneintppd, xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta)
  
  ! Berechnen der D-Matrix
  
  fac = E / (1.d0 - nue**2)
  
  Dmat = 0.d0
  
  k = MAX(1.2d0, 1.d0+0.2d0*area/(25.d0*thickness**2))
  
  Dmat(1,1) =	  fac; Dmat(1,2) = nue*fac;
  Dmat(2,1) = nue*fac; Dmat(2,2) =     fac;
  Dmat(3,3) = E*10.d-6
  Dmat(4,4) = fac * (1.d0 - nue)/2.d0
  Dmat(5,5) = fac * (1.d0 - nue)/(2.d0 * k)
  Dmat(6,6) = fac * (1.d0 - nue)/(2.d0 * k)
  
  !Setzen der H-Matrix
    
  Hmat = 0.d0
  
  Hmat(1,1) = 1.d0
  Hmat(2,5) = 1.d0
  Hmat(3,9) = 1.d0
  Hmat(4,2) = 1.d0
  Hmat(4,4) = 1.d0
  Hmat(5,6) = 1.d0
  Hmat(5,8) = 1.d0
  Hmat(6,3) = 1.d0
  Hmat(6,7) = 1.d0

  quad8_Kgt = 0.d0
  
  do ii=1,numgausspoints
    
    call quad8_ansatzfunction_xieta(xi_vec(ii), eta_vec(ii), Ni=Ni, dNidxi=dNidxi, dNideta=dNideta)
    
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
    
    do jj=1,8
    
      dxdxi   = dxdxi  +dNidxi(jj) * (node_coords(jj,1) + thickness/2.d0*zeta_vec(ii)*normalVec(1,jj))
      dxdeta  = dxdeta +dNideta(jj)* (node_coords(jj,1) + thickness/2.d0*zeta_vec(ii)*normalVec(1,jj))
      dxdzeta = dxdzeta+Ni(jj)*thickness/2.d0*normalVec(1,jj)
    
      dydxi   = dydxi  +dNidxi(jj) * (node_coords(jj,2) + thickness/2.d0*zeta_vec(ii)*normalVec(2,jj))
      dydeta  = dydeta +dNideta(jj)* (node_coords(jj,2) + thickness/2.d0*zeta_vec(ii)*normalVec(2,jj))
      dydzeta = dydzeta+Ni(jj)*thickness/2.d0*normalVec(2,jj)
    
      dzdxi   = dzdxi  +dNidxi(jj) * (node_coords(jj,3) + thickness/2.d0*zeta_vec(ii)*normalVec(3,jj))
      dzdeta  = dzdeta +dNideta(jj)* (node_coords(jj,3) + thickness/2.d0*zeta_vec(ii)*normalVec(3,jj))
      dzdzeta = dzdzeta+Ni(jj)*thickness/2.d0*normalVec(3,jj)
    
    end do
    
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
    
    n2(1:8,1) = 0.5d0 * thickness * dNidxi(1:8)  * zeta_vec(ii)
    n2(1:8,2) = 0.5d0 * thickness * dNideta(1:8) * zeta_vec(ii)
    n2(1:8,3) = 0.5d0 * thickness * Ni(1:8)

    ! Durch MatLab (quad8_iso_egeomat.m) mehr oder weniger automatisch generierter Code, um Operationen mit Null zu vermeiden

    bt( 1: 8, 1) = transMat(1, 1)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt( 9:16, 1) = transMat(2, 1)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt(17:24, 1) = transMat(3, 1)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt(25:32, 1) = transMat(1, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 1) = transMat(1, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 1) = transMat(1, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 1)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 2) = transMat(1, 1)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt( 9:16, 2) = transMat(2, 1)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt(17:24, 2) = transMat(3, 1)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt(25:32, 2) = transMat(1, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 2) = transMat(1, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 2) = transMat(1, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 1)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 3) = transMat(1, 1)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt( 9:16, 3) = transMat(2, 1)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt(17:24, 3) = transMat(3, 1)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt(25:32, 3) = transMat(1, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 3) = transMat(1, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 3) = transMat(1, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 1)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 4) = transMat(1, 2)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt( 9:16, 4) = transMat(2, 2)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt(17:24, 4) = transMat(3, 2)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt(25:32, 4) = transMat(1, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 4) = transMat(1, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 4) = transMat(1, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 2)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 5) = transMat(1, 2)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt( 9:16, 5) = transMat(2, 2)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt(17:24, 5) = transMat(3, 2)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt(25:32, 5) = transMat(1, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 5) = transMat(1, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 5) = transMat(1, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 2)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 6) = transMat(1, 2)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt( 9:16, 6) = transMat(2, 2)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt(17:24, 6) = transMat(3, 2)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt(25:32, 6) = transMat(1, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 6) = transMat(1, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 6) = transMat(1, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 2)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 7) = transMat(1, 3)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt( 9:16, 7) = transMat(2, 3)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt(17:24, 7) = transMat(3, 3)*(jacInv(1, 1)*n1(1:8, 1) + jacInv(1, 2)*n1(1:8, 2))
    bt(25:32, 7) = transMat(1, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 7) = transMat(1, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 7) = transMat(1, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 3)*(jacInv(1, 1)*n2(1:8, 1) + jacInv(1, 2)*n2(1:8, 2) + jacInv(1, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 8) = transMat(1, 3)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt( 9:16, 8) = transMat(2, 3)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt(17:24, 8) = transMat(3, 3)*(jacInv(2, 1)*n1(1:8, 1) + jacInv(2, 2)*n1(1:8, 2))
    bt(25:32, 8) = transMat(1, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 8) = transMat(1, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 8) = transMat(1, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 3)*(jacInv(2, 1)*n2(1:8, 1) + jacInv(2, 2)*n2(1:8, 2) + jacInv(2, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    bt( 1: 8, 9) = transMat(1, 3)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt( 9:16, 9) = transMat(2, 3)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt(17:24, 9) = transMat(3, 3)*(jacInv(3, 1)*n1(1:8, 1) + jacInv(3, 2)*n1(1:8, 2))
    bt(25:32, 9) = transMat(1, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 1) + transMat(2, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 1) + transMat(3, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 1)
    bt(33:40, 9) = transMat(1, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 2) + transMat(2, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 2) + transMat(3, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 2)
    bt(41:48, 9) = transMat(1, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 1, 3) + transMat(2, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 2, 3) + transMat(3, 3)*(jacInv(3, 1)*n2(1:8, 1) + jacInv(3, 2)*n2(1:8, 2) + jacInv(3, 3)*n2(1:8, 3))*nodeTransMat(1:8, 3, 3)
    
    Bmat = MatMul(Hmat, transpose(bt))
    
    Nc = MatMul(Dmat,MatMul(Bmat,reorderedDispl))
    
    lambda(1,1) = Nc(1); lambda(1,2) = Nc(4); lambda(1,3) = Nc(5);
    lambda(2,1) = Nc(4); lambda(2,2) = Nc(2); lambda(2,3) = Nc(6);
    lambda(3,1) = Nc(5); lambda(3,2) = Nc(6); lambda(3,3) = Nc(3);
    
    Sm(1:3,1:3) = lambda
    Sm(4:6,4:6) = lambda
    Sm(7:9,7:9) = lambda
    
    quad8_Kgt = quad8_Kgt + detJac*w_xi(ii)*w_eta(ii)*w_zeta(ii)*MatMul(MatMul(bt,Sm),transpose(bt))
    
  end do
  
  ! sort for nodal degrees of freedom
  ! (u1,u2,u3,u4,u5,u6,u7,u8,v1,v2,v3,v4,...,bez8) -> (u1,v1,w1,bex1,bey1,bez1,...,bez8)
  do ii = 1,48
    do jj = 1,48
      ! uu
      quad8_Kg(index(ii),index(jj)) = quad8_Kgt(ii,jj)
    end do
  end do
  
  !call write_quad_matrix(quad8_Kg, 48,  'Kgel                ', 'beide')
  
  deallocate(xi_vec, eta_vec, zeta_vec)
  deallocate(  w_xi,   w_eta,   w_zeta)
    
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
end subroutine quad8_iso_geostiff
