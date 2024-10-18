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
subroutine quad8_iso_stiff(E,nue,thickness,area,node_coords,quad8_K)

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
  
  double precision, intent(in)			  :: E,nue,thickness,area
  double precision, dimension(8,3), intent(in) 	  :: node_coords
  double precision, dimension(48,48), intent(out) :: quad8_K

  double precision, dimension(8)		  :: dNidxi,dNideta
  double precision                                :: dxdxi, dxdeta, dydxi, dydeta, dzdxi, dzdeta
  double precision, dimension(3,8)      	  :: normalVec 			! Normalenvektoren an den acht Punkten
  double precision, dimension(8,3,3)              :: nodeTransMat		! lokales Koordinatensystem an jedem Knoten
  
  double precision, dimension(8)        	  :: xi_n, eta_n
  
  double precision, dimension(:), allocatable     :: xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta
  
  double precision, dimension(48,48)              :: quad8_Kt			! temporäre Elementsteifigkeitsmatrix
  
  double precision, dimension(3,3,8)              :: nodeTransMat2
  
  double precision, dimension(48,6)               :: bMat
  
  double precision, dimension(3,3)                :: Jac
  
  double precision, dimension(3)                  :: x_loc, y_loc, z_loc	! lokale Koordinatenachsen
  
  double precision                                :: fac, k, detJac
  
  double precision, dimension(6,6)                :: Dmat
  
  double precision                                :: penaltyStiff
  double precision, dimension(8,8)                :: Pmat                                
  double precision, dimension(24,24)              :: extPMat
  double precision, dimension(8,1)                :: vec1
  double precision, dimension(1,8)                :: vec2
  
  integer, parameter                              :: inplaneintppd  = 2		! Integrationspunkte in der Elementebene
  integer, parameter                              :: outplaneintppd = 2		! Integrationspunkte in Dickenrichtung
  
  integer                                         :: ii, jj
  integer                                         :: numgausspoints
  double precision, dimension(3,3)                :: transMat
  
  integer, dimension(48), parameter               :: index = (/ 1, 7,13,19,25,31,37,43, 2, 8,14,20,26,32,38,44, 3, 9,15,21,27,33,39,45, 4,10,16,22,28,34,40,46, 5,11,17,23,29,35,41,47, 6,12,18,24,30,36,42,48/)

  xi_n(1) = -1.d0; eta_n(1) = -1.d0
  xi_n(2) =  1.d0; eta_n(2) = -1.d0
  xi_n(3) =  1.d0; eta_n(3) =  1.d0
  xi_n(4) = -1.d0; eta_n(4) =  1.d0
  xi_n(5) =  0.d0; eta_n(5) = -1.d0
  xi_n(6) =  1.d0; eta_n(6) =  0.d0
  xi_n(7) =  0.d0; eta_n(7) =  1.d0
  xi_n(8) = -1.d0; eta_n(8) =  0.d0
  
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
  
  Dmat(1,1) = fac;       Dmat(1,2) = nue*fac;
  Dmat(2,1) = Dmat(1,2); Dmat(2,2) =     fac;
  Dmat(3,3) = E*10.d-6
  Dmat(4,4) = fac * (1.d0 - nue)/2.d0
  Dmat(5,5) = fac * (1.d0 - nue)/(2.d0 * k)
  Dmat(6,6) = Dmat(5,5)

  quad8_Kt = 0.d0
  
  do ii=1,numgausspoints
    
    call quad8_bmat(xi_vec(ii), eta_vec(ii), zeta_vec(ii), node_coords, thickness, normalVec, nodeTransMat, bmat, detJac, Jac, transMat)

    quad8_Kt = quad8_Kt + detJac*w_xi(ii)*w_eta(ii)*w_zeta(ii)*MatMul(MatMul(Bmat,Dmat),transpose(Bmat))

  end do

  penaltyStiff = 1.d0/(1.d0-nue**2)/100000.d0*thickness*area*E
  
  pMat = -penaltyStiff/7.d0
  
  do ii = 1,8
    pMat(ii,ii) = penaltyStiff
  end do
  
  do ii = 1,3
    vec1(1:8,1) = nodeTransMat2(ii,3,1:8)
    do jj = 1,3
      vec2(1,1:8) = nodeTransMat2(jj,3,1:8)
      extPMat((ii-1)*8+1:ii*8,(jj-1)*8+1:jj*8) = matmul(vec1,vec2)*pMat
    end do
  end do
    
  quad8_Kt(25:48,25:48) = quad8_Kt(25:48,25:48) + extPMat
  
  ! sort for nodal degrees of freedom
  ! (u1,u2,u3,u4,u5,u6,u7,u8,v1,v2,v3,v4,...,bez8) -> (u1,v1,w1,bex1,bey1,bez1,...,bez8)
  do ii = 1,48
    do jj = 1,48
      ! uu
      quad8_K(index(ii),index(jj)) = quad8_Kt(ii,jj)
    end do
  end do

  deallocate(xi_vec, eta_vec, zeta_vec)
  deallocate(  w_xi,   w_eta,   w_zeta)

end subroutine quad8_iso_stiff
