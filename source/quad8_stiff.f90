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
!> Calculates element stiffness matrix for layered 8-node shell element quad8
!
!> @details
!> Calculation of element stiffness matrix for layered 8-node shell element
!> quad8; integration is done layerwise, zeta = [-1,1] is substituted for
!> each layer to run over thickness of layer such that the same integration point
!> values can be used in every layer; by this, the layer thicknesses are taken into
!> account as relative to the total laminat thickness which therefore needs not to be
!> equal to the actual thickness of the element as layer thicknesses are scaled to it;
!
!> @author
!> Andreas Hauffe, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 23.05.2017
!
!> $Id: quad8_stiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_stiff(nlay,Clay,Tlay,lth,tth,pStiff,node_coords,nip,nop,opit,area,quad8Num,elemNum,quad8_K)
! =================================================================================================
! use
!
  use integration_schemes
  
  implicit none
!
! Input
!
  integer, intent(in)                               :: nlay               !< Total number of layers in current element
  double precision, intent(in)                      :: Clay(6,6,nlay)     !< (6,6,nlay)-Array which contains the material stiffness matrix in the element coordinate system for each layer
  double precision, intent(in)                      :: Tlay(6,6,nlay)     !< (6,6,nlay)-Array which contains the Transformation matrix in the element coordinate system for each layer
  double precision, intent(in)                      :: lth(nlay)          !< (nlay)-Array which contains the thickness for each layer
  double precision, intent(in)                      :: tth                !< Total thickness of laminat (sum of all values in lth)
  double precision, intent(in)                      :: pStiff             !< (ym11+ym22+ym33)/3.d0
  double precision, intent(in)                      :: node_coords(8,3)   !< (8,3)-Array containing global coordinates of element nodes
  integer, intent(in)                               :: nip                !< Number of integration points in plane per direction
  integer, intent(in)                               :: nop(nlay)          !< (nlay)-Array with number of integration points over layer thickness for each layer
  integer, intent(in)                               :: opit               !< Out-Off-Plane integration type
  double precision, intent(in)                      :: area               !< element reference surface area
  integer, intent(in)                               :: quad8Num           !< number of the quad8 element (only used for error messages)
  integer, intent(in)                               :: elemNum            !< element number of the quad8 element (only used for error messages)
!
! Output
!
  double precision, dimension(48,48), intent(out)   :: quad8_K            !< 48x48-Element stiffness matrix (8x6 dof)
!
! Internal
!
  double precision, dimension(3,8)                  :: normalVec                ! Normalenvektoren an den acht Punkten
  double precision, dimension(8,3,3)                :: nodeTransMat             ! lokales Koordinatensystem an jedem Knoten
  double precision, dimension(:), allocatable       :: xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta
  double precision, dimension(48,48)                :: quad8_Kt                 ! temporäre Elementsteifigkeitsmatrix
  double precision, dimension(3,3,8)                :: nodeTransMat2
  double precision, dimension(48,6)                 :: bMat
  double precision, dimension(6,6)                  :: Cloc
  double precision, dimension(3,3)                  :: Jac
  double precision, dimension(3)                    :: xElem
  double precision                                  :: k, detJac
  double precision                                  :: penaltyStiff
  double precision, dimension(8,8)                  :: Pmat                                
  double precision, dimension(24,24)                :: extPMat
  double precision, dimension(8,1)                  :: vec1
  double precision, dimension(1,8)                  :: vec2
  integer, parameter                                :: outplaneintppd = 2       ! Integrationspunkte in Dickenrichtung
  integer                                           :: ii, jj
  integer                                           :: numIntPointsIP, numIntPointsOoP
  integer                                           :: iLay
  double precision, dimension(3,3)                  :: transMat
  double precision                                  :: zeta_elem, lth_sum
  
  integer                                           :: err_code=0
  
  ! Die eps_zRotStiff (0.99999d0) sind per Reverse Engineering entstanden, um die selbe Steifigkeitsmatrix wie das Shell91 zu erhalten
  double precision, parameter                       :: eps_zRotStiff = 0.99999d0
  
  integer, dimension(48), parameter                 :: index = (/ 1, 7,13,19,25,31,37,43, 2, 8,14,20,26,32,38,44, 3, 9,15,21,27,33,39,45, 4,10,16,22,28,34,40,46, 5,11,17,23,29,35,41,47, 6,12,18,24,30,36,42,48/)
  
  ! Berechnung der Normalen auf der Elementmittelfläche und der dazugehörigen Koordinatensysteme mit Transformationsmatrizen
  call quad8_normals(node_coords,normalVec,nodeTransMat,nodeTransMat2)

  !*****************************************************
  ! Integration
  !*****************************************************
  
  quad8_Kt = 0.d0
  
  xElem(1) = node_coords(2,1) - node_coords(1,1)
  xElem(2) = node_coords(2,2) - node_coords(1,2)
  xElem(3) = node_coords(2,3) - node_coords(1,3)
  
  numIntPointsIP = nip*nip
    
  allocate(xi_vec(numIntPointsIP), eta_vec(numIntPointsIP))
  allocate(  w_xi(numIntPointsIP),   w_eta(numIntPointsIP))
  
  call integration_points_2d(nip, "Gauss",  xi_vec, w_xi, &
                             nip, "Gauss", eta_vec, w_eta)
  
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
            
            ! Beschleunigte Berechung der B-Matrix durch Verzicht auf viele der Nulloperationen
            call quad8_bmat(xi_vec(ii), eta_vec(ii), zeta_elem, node_coords, tth, normalVec, nodeTransMat, bmat, detJac, Jac, transMat)

            quad8_Kt = quad8_Kt + detJac*w_xi(ii)*w_eta(ii)*w_zeta(jj)*lth(iLay)/tth*MatMul(MatMul(Bmat,Cloc),transpose(Bmat))
            
            if (detJac < 0.d0) then
              write(*,*) 'for quad8 number:                     ', quad8Num
              write(*,*) 'element number:                       ', elemNum
              write(*,*) 'negative determinant of the jacobian: ', detJac
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

  penaltyStiff = 0.d0
  do iLay = 1,nlay 
    penaltyStiff = penaltyStiff + pStiff*lth(iLay)*area/100000.d0
  end do
  
  pMat = -penaltyStiff/7.d0*eps_zRotStiff
  
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

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_stiff'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_stiff
