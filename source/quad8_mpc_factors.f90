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
!
!> @details
!
!> @author 
!
!> $Id: quad8_mpc_factors.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_mpc_factors(node_coords, xi, eta, factors)

  use mat_func

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
  
  double precision, intent(in)                    :: xi, eta
  double precision, dimension(8,3), intent(in)    :: node_coords
  double precision, dimension(3,24), intent(out)  :: factors

  double precision, dimension(8)                  :: Ni, dNidxi,dNideta
  double precision                                :: dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta, dzdxi, dzdeta, dzdzeta
  double precision, dimension(3,8)                :: normalVec              ! Normalenvektoren an den acht Punkten
  double precision, dimension(3)                  :: x_loc, y_loc, z_loc    ! lokale Koordinatenachsen
  
  double precision, dimension(8)                  :: xi_n, eta_n
  
    
  double precision, dimension(3,3)                :: jac, jacInv            ! Jakobi-Matrix und inverse Jakobi-Matrix an den Integrationspunkten
  
  double precision, dimension(3,3)                :: transMat
    
  double precision, dimension(3,8)                :: n1
  
  double precision, dimension(9,24)               :: bt
  
  double precision, dimension(3,24)               :: BMat
  
  double precision, dimension(9,9)                :: expTransMat
  
  double precision, dimension(3,9)                :: Hmat
    
  double precision                                :: detJac
    
  double precision                                :: zeta, thickness

  integer                                         :: ii, jj, mm, nn, kk
    
  logical                                             :: ok_flag

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
    
  end do
  
  !Setzen der H-Matrix
    
  Hmat = 0.d0
  
  Hmat(1,8) =  1.d0
  Hmat(2,7) = -1.d0
  Hmat(3,2) = -.5d0
  Hmat(3,4) =  .5d0

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
  
  thickness = 1.d0
  zeta = 0.d0
  
  do jj=1,8
  
    dxdxi   = dxdxi  +dNidxi(jj) * (node_coords(jj,1) + thickness/2.d0*zeta*normalVec(1,jj))
    dxdeta  = dxdeta +dNideta(jj)* (node_coords(jj,1) + thickness/2.d0*zeta*normalVec(1,jj))
    dxdzeta = dxdzeta+Ni(jj)*thickness/2.d0*normalVec(1,jj)
  
    dydxi   = dydxi  +dNidxi(jj) * (node_coords(jj,2) + thickness/2.d0*zeta*normalVec(2,jj))
    dydeta  = dydeta +dNideta(jj)* (node_coords(jj,2) + thickness/2.d0*zeta*normalVec(2,jj))
    dydzeta = dydzeta+Ni(jj)*thickness/2.d0*normalVec(2,jj)
  
    dzdxi   = dzdxi  +dNidxi(jj) * (node_coords(jj,3) + thickness/2.d0*zeta*normalVec(3,jj))
    dzdeta  = dzdeta +dNideta(jj)* (node_coords(jj,3) + thickness/2.d0*zeta*normalVec(3,jj))
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
  
  do mm = 1,3
    transMat(mm,1) = x_loc(mm)
    transMat(mm,2) = y_loc(mm)
    transMat(mm,3) = z_loc(mm)
  end do
  
  do jj = 1,8
  
    n1(1,jj) = dNidxi(jj)
    n1(2,jj) = dNideta(jj)
    n1(3,jj) = 0.d0
  
  end do
  
  jacInv = matmul(transpose(transMat), jacInv)
  
  n1 = matmul(jacInv, n1)  ! T^T * J^-1 * n_1

  bt = 0.d0

  bt(1:3, 1: 8) = n1
  bt(4:6, 9:16) = n1
  bt(7:9,17:24) = n1
  
  do mm = 1,3
    do nn = 1,3
      do kk = 1,3
        expTransMat((nn-1)*3+kk,(mm-1)*3+kk) = transMat(mm,nn)
      end do
    end do
  end do
  
  bt = MatMul(expTransMat, bt)
  
  Bmat = MatMul(Hmat, bt)
  
  factors = MatMul(transMat,Bmat)
  
end subroutine quad8_mpc_factors
