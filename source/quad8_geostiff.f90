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
!> Calculates element geometric stiffness matrix for layered 8-node shell element quad8
!
!> @details
!> Calculation of element geometric stiffness matrix for layered 8-node shell element
!> quad8; integration is done layerwise, zeta = [-1,1] is substituted for
!> each layer to run over thickness of layer such that the same integration point
!> values can be used in every layer;
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
!> $Id: quad8_geostiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
subroutine quad8_geostiff(nlay,Clay,Tlay,lth,tth,node_coords,nip,nop,opit,area,quad8_Kg,disp_q8,ntemps,ath,is_temperature,is_quad8temp)
! =================================================================================================
! use
!
  use integration_schemes
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
!
! Input
!
  integer, intent(in)                               :: nlay             !< Total number of layers in current element
  double precision, intent(in)                      :: Clay(6,6,nlay)   !< (6,6,nlay)-Array which contains the material stiffness matrix in the element coordinate system for each layer
  double precision, intent(in)                      :: Tlay(6,6,nlay)   !< (6,6,nlay)-Array which contains the Transformation matrix in the element coordinate system for each layer
  double precision, intent(in)                      :: lth(nlay)        !< (nlay)-Array which contains the thickness for each layer
  double precision, intent(in)                      :: tth              !< Total thickness of laminat (sum of all values in lth)
  double precision, intent(in)                      :: node_coords(8,3) !< (8,3)-Array containing global coordinates of element nodes
  integer, intent(in)                               :: nip              !< Number of integration points in plane per direction
  integer, intent(in)                               :: nop(nlay)        !< (nlay)-Array with number of integration points over layer thickness for each layer
  integer, intent(in)                               :: opit             !< Out-Off-Plane integration type
  double precision, intent(in)                      :: area             !< element reference surface area
  double precision, dimension(48), intent(in)       :: disp_q8          !< global displacements for each node
  double precision, intent(in)                      :: ntemps(8)        !< (8)-Array containing the temperature at each node
  double precision, intent(in)                      :: ath(6,nlay)      !< (6,nlay)-Array containing the coefficients of thermal expansion
  logical, intent(in)                               :: is_temperature
  logical, intent(in)                               :: is_quad8temp
!
! Output
!
  double precision, dimension(48,48), intent(out)   :: quad8_Kg         !< 48x48-Element geometric stiffness matrix (8x6 dof)

  double precision, dimension(8)                    :: Ni, dNidxi,dNideta
  double precision                                  :: dxdxi, dxdeta, dxdzeta, dydxi, dydeta, dydzeta, dzdxi, dzdeta, dzdzeta
  double precision, dimension(3,8)                  :: normalVec                ! Normalenvektoren an den acht Punkten
  double precision, dimension(8,3,3)                :: nodeTransMat             ! lokales Koordinatensystem an jedem Knoten
  
  double precision, dimension(:), allocatable       :: xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta
  
  double precision, dimension(48,48)                :: quad8_Kgt                ! temporäre Elementsteifigkeitsmatrix
  
  double precision, dimension(48)                   :: reorderedDispl
  
  double precision, dimension(6)                    :: Nc
  
  double precision, dimension(3,3)                  :: lambda
  
  double precision, dimension(9,9)                  :: Sm
  
  double precision, dimension(3,3)                  :: jac, jacInv              ! Jakobi-Matrix und inverse Jakobi-Matrix an den Integrationspunkten
  double precision, dimension(3)                    :: xElem

  double precision, dimension(3,3)                  :: transMat
  
  double precision, dimension(3,3,8)                :: nodeTransMat2
  
  double precision, dimension(8,3)                  :: n1, n2
  
  double precision, dimension(48,9)                 :: bt
  
  double precision, dimension(6,9)                  :: Hmat
  
  double precision, dimension(6,48)                 :: bMat
  
  double precision, dimension(6,6)                  :: Cloc
  
  double precision, dimension(6)                    :: totstrain, eps_th, mechstrain
  
  double precision                                  :: k, detJac
  integer                                           :: numIntPointsIP, numIntPointsOoP
  double precision                                  :: zeta_elem, lth_sum
  integer                                           :: iLay
  
  integer                                           :: ii, jj
  
  integer, dimension(48), parameter                 :: index = (/ 1, 7,13,19,25,31,37,43, 2, 8,14,20,26,32,38,44, 3, 9,15,21,27,33,39,45, 4,10,16,22,28,34,40,46, 5,11,17,23,29,35,41,47, 6,12,18,24,30,36,42,48/)

  logical                                           :: ok_flag
  
  integer                                           :: err_code=0
  
  do jj = 1,6
    do ii = 1,8
      reorderedDispl((jj-1)*8+ii) = disp_q8((ii-1)*6+jj)
    end do
  end do
  
  ! Berechnung der Normalen auf der Elementmittelfläche und der dazugehörigen Koordinatensysteme mit Transformationsmatrizen
  call quad8_normals(node_coords,normalVec,nodeTransMat,nodeTransMat2)

  !*****************************************************
  ! Integration
  !*****************************************************
  
  quad8_Kgt = 0.d0
  
  xElem(1) = node_coords(2,1) - node_coords(1,1)
  xElem(2) = node_coords(2,2) - node_coords(1,2)
  xElem(3) = node_coords(2,3) - node_coords(1,3)
  
  numIntPointsIP = nip*nip
    
  allocate(xi_vec(numIntPointsIP), eta_vec(numIntPointsIP))
  allocate(  w_xi(numIntPointsIP),   w_eta(numIntPointsIP))
  
  call integration_points_2d(nip, "Gauss",  xi_vec, w_xi, &
                             nip, "Gauss", eta_vec, w_eta)
  
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
  
  do ii = 1,numIntPointsIP
  
    ! Berechnung des Koordinatensystems im Integrationspunkt (nur in der Mittelebene, gilt für die gesamte Dickenrichtung)
    call quad8_ip_coordsys(node_coords, xi_vec(ii), eta_vec(ii), xElem, transMat)
  
    lth_sum = 0.d0
  
    ! --- Beginn der Schleife über die Lagen

    do iLay = 1,nlay 
        
        ! Integrations in Dickenrichtung
        ! Shell 91 - 3-Punkt-Simson
        ! Shell 93 - 2-Punkt-Gauss

        numIntPointsOoP = nop(iLay)
        
        allocate(zeta_vec(numIntPointsOoP))
        allocate(  w_zeta(numIntPointsOoP))
        
        if (opit .EQ. 93) then
          call gauss_integration(numIntPointsOoP, zeta_vec, w_zeta)
        else if (opit .EQ. 91) then
          call simpson_integration(numIntPointsOoP, zeta_vec, w_zeta)
        else
          write(*,*) 'Unknown integration type'
          err_code = 2
          goto 9999
        end if

        ! calculate summed thickness until current layer
        lth_sum = lth_sum + lth(iLay)

        ! Berechnen der C-Matrix
        k = MAX(1.2d0, 1.d0+0.2d0*area/(25.d0*tth**2))
        
        Cloc(1:6,1:6) = Clay(1:6,1:6,iLay)
        
        Cloc(5,5) = Cloc(5,5) / k
        Cloc(6,6) = Cloc(6,6) / k
        
        Cloc = matmul(transpose(Tlay(1:6,1:6,iLay)),matmul(Cloc,Tlay(1:6,1:6,iLay)))
        
        do jj=1,numIntPointsOoP

            ! transform layerwise zeta-value of integration point zeta_ip(ii)
            ! to elementwise zeta_value zeta_elem
            ! (see eq.(22) in Panda, S.; Natarajan, R.: Analysis of laminated composite shell
            !  structures by finite element method. Computers & Structures 14(3-4). 1981. pp. 225-230)
            !zeta_elem = -1.d0 + (2.d0*lth_sum - lth(iLay)*(1.d0-zeta_vec(ii)))/tth
            ! angepasst, da obere Variante durch Rundung unsymmetisch bezüglich null werden kann
            zeta_elem = (2.d0*lth_sum - lth(iLay) - tth)/tth + lth(iLay)/tth*zeta_vec(jj)
            
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
            
            dxdxi   = SUM(dNidxi(1:8) * (node_coords(1:8,1) + tth/2.d0*zeta_elem*normalVec(1,1:8)))
            dxdeta  = SUM(dNideta(1:8)* (node_coords(1:8,1) + tth/2.d0*zeta_elem*normalVec(1,1:8)))
            dxdzeta = SUM(Ni(1:8)*tth/2.d0*normalVec(1,1:8))
            
            dydxi   = SUM(dNidxi(1:8) * (node_coords(1:8,2) + tth/2.d0*zeta_elem*normalVec(2,1:8)))
            dydeta  = SUM(dNideta(1:8)* (node_coords(1:8,2) + tth/2.d0*zeta_elem*normalVec(2,1:8)))
            dydzeta = SUM(Ni(1:8)*tth/2.d0*normalVec(2,1:8))
            
            dzdxi   = SUM(dNidxi(1:8) * (node_coords(1:8,3) + tth/2.d0*zeta_elem*normalVec(3,1:8)))
            dzdeta  = SUM(dNideta(1:8)* (node_coords(1:8,3) + tth/2.d0*zeta_elem*normalVec(3,1:8)))
            dzdzeta = SUM(Ni(1:8)*tth/2.d0*normalVec(3,1:8))
            
            jac(1,1) =   dxdxi; jac(1,2) =   dydxi; jac(1,3) =   dzdxi;
            jac(2,1) =  dxdeta; jac(2,2) =  dydeta; jac(2,3) =  dzdeta;
            jac(3,1) = dxdzeta; jac(3,2) = dydzeta; jac(3,3) = dzdzeta;
            
            call m33inv(jac, jacInv, detJac, ok_flag)
            
            jacInv = matmul(transpose(transMat), jacInv)
            
            n1(1:8,1) = dNidxi(1:8)
            n1(1:8,2) = dNideta(1:8)
            n1(1:8,3) = 0.d0
            
            n2(1:8,1) = 0.5d0 * tth * dNidxi(1:8)  * zeta_elem
            n2(1:8,2) = 0.5d0 * tth * dNideta(1:8) * zeta_elem
            n2(1:8,3) = 0.5d0 * tth * Ni(1:8)

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
            
            totstrain(1:6) = MatMul(Bmat,reorderedDispl)
            
            if ((is_temperature == .true.) .or. (is_quad8temp == .true.))then
                eps_th(1:6)  = SUM(ntemps(1:8)*Ni(1:8)) * ath(1:6,iLay)
                mechstrain(1:6) = totstrain(1:6) - eps_th(1:6)
            else
                mechstrain = totstrain
            end if
            mechstrain(3) = 0.d0

            Nc = MatMul(Cloc,mechstrain)
            
            lambda(1,1) = Nc(1); lambda(1,2) = Nc(4); lambda(1,3) = Nc(5);
            lambda(2,1) = Nc(4); lambda(2,2) = Nc(2); lambda(2,3) = Nc(6);
            lambda(3,1) = Nc(5); lambda(3,2) = Nc(6); lambda(3,3) = Nc(3);
            
            Sm(1:3,1:3) = lambda
            Sm(4:6,4:6) = lambda
            Sm(7:9,7:9) = lambda
            
            quad8_Kgt = quad8_Kgt + detJac*w_xi(ii)*w_eta(ii)*w_zeta(jj)*lth(iLay)/tth*MatMul(MatMul(bt,Sm),transpose(bt))
            
            if (detJac < 0.d0) then
              err_code = 2
              goto 9999
            end if

        end do
        
        deallocate(zeta_vec)
        deallocate(  w_zeta)
    
    end do

    ! --- Ende der Schleife über die Lagen
  
  end do

  deallocate(xi_vec, eta_vec)
  deallocate(  w_xi,   w_eta)
  
  ! sort for nodal degrees of freedom
  ! (u1,u2,u3,u4,u5,u6,u7,u8,v1,v2,v3,v4,...,bez8) -> (u1,v1,w1,bex1,bey1,bez1,...,bez8)
  do ii = 1,48
    do jj = 1,48
      ! uu
      quad8_Kg(index(ii),index(jj)) = quad8_Kgt(ii,jj)
    end do
  end do
  
  !call write_quad_matrix(quad8_Kg, 48,  'Kgel                ', 'beide')

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_geostiff'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_geostiff
