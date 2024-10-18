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
!> Calculation of the B-Matrix for an 8-node quadratic shell element
!
!> @details
!> 
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 29.06.2013 
!
!> $Id: quad8_bmat.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_bmat(xi, eta, zeta, node_coords, thickness, normalVec, nodeTransMat, bmat, detJac, Jac, transMat)
! =================================================================================================
! use
!
  use mat_func

  implicit none
!
! =================================================================================================
!
! Interface
!
  interface 
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
      double precision, intent(in)                                  :: xi, eta
      
      double precision, dimension(8), optional, intent(out)         :: Ni
      double precision, dimension(8), optional, intent(out)         :: dNidxi, dNideta
      double precision, dimension(8), optional, intent(out)         :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
  end interface
  
    
  double precision, intent(in)                      :: xi               !< xi coordinate inside the element (-1 <= xi <= 1)
  double precision, intent(in)                      :: eta              !< eta coordinate inside the element (-1 <= eta <= 1)
  double precision, intent(in)                      :: zeta             !< zeta coordinate inside the element (-1 <= zeta <= 1)
  double precision, dimension(8,3), intent(in)      :: node_coords      !< (8,3)-Array containing global coordinates of element nodes
  double precision, intent(in)                      :: thickness        !< total thickness of the shell element
  double precision, dimension(8,3,3), intent(in)    :: nodeTransMat     !< (8,3,3)-Array containing the transformation matrix of each node for calculation part of the displacements out of the rotations
  
  double precision, dimension(48,6), intent(out)    :: bMat             !< (48,6)-Array containing the B matrix which transforms the global nodal displacements to the strain vector
  double precision, dimension(3,3), intent(out)     :: jac              !< (3,3)-Array containing the Jacobi matrix at the integration point
  double precision, intent(out)                     :: detJac           !< determinant of the jacobi matrix

  double precision, dimension(8)                    :: Ni, dNidxi,dNideta
  double precision                                  :: dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta, dzdxi, dzdeta, dzdzeta
  double precision, dimension(3,8)                  :: normalVec        ! Normalenvektoren an den acht Punkten
  double precision, dimension(3,3)                  :: jacInv           ! inverse Jakobi-Matrix am Integrationspunkt
  double precision, dimension(3,3), intent(in)      :: transMat
  double precision, dimension(8,3)                  :: n1, n2

  logical                                               :: ok_flag

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
    
end subroutine quad8_bmat
