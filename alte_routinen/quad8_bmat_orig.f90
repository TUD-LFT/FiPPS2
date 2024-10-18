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
subroutine quad8_bmat_orig(xi, eta, zeta, node_coords, thickness, normalVec, nodeTransMat, bmatret, detJac, Jac, transMat)

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
  
    
  double precision, intent(in)                      :: xi, eta, zeta
  double precision, dimension(8,3), intent(in)      :: node_coords
  double precision, intent(in)                      :: thickness
  double precision, dimension(8,3,3), intent(in)    :: nodeTransMat     ! lokales Koordinatensystem an jedem Knoten
  double precision, dimension(3,8), intent(in)      :: normalVec        ! Normalenvektoren an den acht Punkten
  
  double precision, dimension(48,6), intent(out)    :: bMatRet
  double precision, dimension(3,3)                  :: jac, jacInv              ! Jakobi-Matrix am Integrationspunkt
  double precision, intent(out)                     :: detJac
  
  double precision, dimension(6,9)                  :: Hmat
  double precision, dimension(8)                    :: Ni, dNidxi, dNideta
  double precision                                  :: dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta, dzdxi, dzdeta, dzdzeta
  double precision, dimension(3)                    :: x_loc, y_loc, z_loc  ! lokale Koordinatenachsen
  double precision, dimension(3,3), intent(out)     :: transMat
  double precision, dimension(3,8)                  :: n1, n2
  double precision, dimension(9,48)                 :: bt
  double precision, dimension(6,48)                 :: bMat
  
  double precision, dimension(9,24)                 :: n1Exp, n2Exp
  double precision, dimension(9,9)                  :: JacInvExp, transMatExp, transMatExp2
  
  integer                                           :: ii, jj
  
  logical                                           :: ok_flag
  
  double precision, dimension(24,24)                :: t1 ! Tranformation der globalen Winkel der Knoten in das lokale System
  
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
  
  t1 = 0.d0
  
  ! Transformation der globalen Rotationen ins Knotenkoordinatensystem
  do jj = 1,8
  
    t1(jj   ,jj   ) = nodeTransMat(jj,1,1)
    t1(jj   ,jj+ 8) = nodeTransMat(jj,1,2)
    t1(jj   ,jj+16) = nodeTransMat(jj,1,3)
    t1(jj+ 8,jj   ) = nodeTransMat(jj,2,1)
    t1(jj+ 8,jj+ 8) = nodeTransMat(jj,2,2)
    t1(jj+ 8,jj+16) = nodeTransMat(jj,2,3)
    t1(jj+16,jj   ) = nodeTransMat(jj,3,1)
    t1(jj+16,jj+ 8) = nodeTransMat(jj,3,2)
    t1(jj+16,jj+16) = nodeTransMat(jj,3,3)
  
  end do
  
  call quad8_ansatzfunction_xieta(xi, eta, Ni=Ni, dNidxi=dNidxi, dNideta=dNideta)
  
  do jj = 1,8
  
    n1(1,jj) = dNidxi(jj)
    n1(2,jj) = dNideta(jj)
    n1(3,jj) = 0.d0
    
    n2(1,jj) = 0.5d0 * thickness * dNidxi(jj)  * zeta
    n2(2,jj) = 0.5d0 * thickness * dNideta(jj) * zeta
    n2(3,jj) = 0.5d0 * thickness * Ni(jj)
  
  end do
  
  n1Exp = 0.d0
  n2Exp = 0.d0
  bt = 0.d0
  
  ! Glg. 1.12
  n1Exp(1:3, 1: 8) = n1(1:3,1:8)
  n1Exp(4:6, 9:16) = n1(1:3,1:8)
  n1Exp(7:9,17:24) = n1(1:3,1:8)
  
  n2Exp(1:3, 1: 8) = n2(1:3,1:8)
  n2Exp(4:6, 9:16) = n2(1:3,1:8)
  n2Exp(7:9,17:24) = n2(1:3,1:8)
  
  n2Exp = matmul(n2Exp,t1)
  
  ! Glg. 1.13
  bt(1:9, 1:24) = n1Exp(1:9,1:24)
  bt(1:9,25:48) = n2Exp(1:9,1:24)
  
  ! Aufstellen der Jakobimatrix
  
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
  
  JacInvExp = 0.d0
  
  ! Glg. 1.14
  JacInvExp(1:3,1:3) = jacInv(1:3,1:3)
  JacInvExp(4:6,4:6) = jacInv(1:3,1:3)
  JacInvExp(7:9,7:9) = jacInv(1:3,1:3)
  
  bt = matmul(JacInvExp,bt)

  
  ! Aufstellen der Transformationsmatrix für das Gauspunktkoordinatensystem
  dxdxi  = SUM(dNidxi(1:8) *node_coords(1:8,1))
  dxdeta = SUM(dNideta(1:8)*node_coords(1:8,1))
  
  dydxi  = SUM(dNidxi(1:8) *node_coords(1:8,2))
  dydeta = SUM(dNideta(1:8)*node_coords(1:8,2))
  
  dzdxi  = SUM(dNidxi(1:8) *node_coords(1:8,3))
  dzdeta = SUM(dNideta(1:8)*node_coords(1:8,3))
  
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
  y_loc = y_loc / SQRT(DOT_PRODUCT(y_loc,y_loc))
   
  ! cross_prod(z_vec, x_vec) (sollte normiert sein, da x und z normiert wurden)
  x_loc(1) = y_loc(2)*z_loc(3) - z_loc(2)*y_loc(3)
  x_loc(2) = y_loc(3)*z_loc(1) - z_loc(3)*y_loc(1)
  x_loc(3) = y_loc(1)*z_loc(2) - z_loc(1)*y_loc(2)
  x_loc = x_loc / SQRT(DOT_PRODUCT(x_loc,x_loc))
  
  transMat(1:3,1) = x_loc(1:3)
  transMat(1:3,2) = y_loc(1:3)
  transMat(1:3,3) = z_loc(1:3)
  
  !
  transMatExp = 0.d0
  transMatExp(1:3,1:3) = transMat(1:3,1:3)
  transMatExp(4:6,4:6) = transMat(1:3,1:3)
  transMatExp(7:9,7:9) = transMat(1:3,1:3)
  
  transMat = transpose(transMat)
  
  transMatExp2 = 0.d0
  do ii = 1,3
    transMatExp2(ii  ,ii  ) = transMat(1,1)
    transMatExp2(ii  ,ii+3) = transMat(1,2)
    transMatExp2(ii  ,ii+6) = transMat(1,3)
    transMatExp2(ii+3,ii  ) = transMat(2,1)
    transMatExp2(ii+3,ii+3) = transMat(2,2)
    transMatExp2(ii+3,ii+6) = transMat(2,3)
    transMatExp2(ii+6,ii  ) = transMat(3,1)
    transMatExp2(ii+6,ii+3) = transMat(3,2)
    transMatExp2(ii+6,ii+6) = transMat(3,3)
  end do    
  
  ! Glg. 1.17
  bt = matmul(transMatExp2,matmul(transpose(transMatExp),bt))
  
  ! Glg. 1.21
  Bmat = MatMul(Hmat, bt)
  
  BmatRet = transpose(bmat)
    
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
      
end subroutine quad8_bmat_orig
