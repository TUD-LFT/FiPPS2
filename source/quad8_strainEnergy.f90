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
!> Calculates elementwise total strain energy for quad8 elements
!
!> @details
!> Elementwise calculation of the total strain energy (tse) for quad8 elements. The total
!> strain energy is calculated as the strain tensor multiplied with the stress tensor, multiplied
!> by 0.5 and integrated over the element volume. Integration over the element volume is performed
!> numerically by using the Gauss (non-layered option) or Simpson (layered option) rule out of
!> plane and Gauss rule in plane.
!
!> @author
!> Florian Dexl, TU Dresden, wiss. Mitarbeiter
!
!> @date
!> 01.07.2020
!
!> $Id: quad8s_temploadvec.f90 293 2018-02-07 12:51:39Z DOM\ahauffe $
!> $Author: DOM\ahauffe $
!> $Revision: 293 $
!> $Date: 2018-02-07 13:51:39 +0100 (Mi, 07. Feb 2018) $
!
! =================================================================================================
subroutine quad8_strainEnergy(nlay,Clay,Tlay,lth,tth,node_coords,nip,nop,opit,area,quad8_tse,disp_q8,ntemps,ath,is_temperature,is_quad8temp)
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
  double precision, intent(out)                     :: quad8_tse        !< Total strain energy of quad8 element
  double precision, dimension(3,8)                  :: normalVec                ! Normalenvektoren an den acht Punkten
  double precision, dimension(8,3,3)                :: nodeTransMat             ! lokales Koordinatensystem an jedem Knoten
  double precision, dimension(:), allocatable       :: xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta
  double precision, dimension(48)                   :: reorderedDispl
  double precision, dimension(6)                    :: Nc
  double precision, dimension(8)                    :: Ni
  double precision, dimension(3,3)                  :: jac
  double precision, dimension(3)                    :: xElem
  double precision, dimension(3,3)                  :: transMat
  double precision, dimension(3,3,8)                :: nodeTransMat2
  double precision, dimension(48,6)                 :: bMat
  double precision, dimension(6,6)                  :: Cloc
  double precision, dimension(6)                    :: totstrain, mechstrain  
  double precision                                  :: k, detJac
  integer                                           :: numIntPointsIP, numIntPointsOoP
  double precision                                  :: zeta_elem, lth_sum
  integer                                           :: iLay
  integer                                           :: ii, jj
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
  
  quad8_tse = 0.d0
  
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
    
            call quad8_ansatzfunction_xieta(xi_vec(ii), eta_vec(ii), Ni=Ni)
            
            call quad8_bmat(xi_vec(ii), eta_vec(ii), zeta_elem, node_coords, tth, normalVec, nodeTransMat, bmat, detJac, Jac, transMat)
            
            totstrain(1:6) = MatMul(transpose(Bmat),reorderedDispl)
            
            if ((is_temperature == .true.) .or. (is_quad8temp == .true.))then
                mechstrain(1:6) = totstrain(1:6) - SUM(ntemps(1:8)*Ni(1:8)) * ath(1:6,iLay)
            else
                mechstrain = totstrain
            end if
            mechstrain(3) = 0.d0
            
            Nc = MatMul(Cloc,totstrain)
            
            quad8_tse = quad8_tse + detJac*w_xi(ii)*w_eta(ii)*w_zeta(jj)*lth(iLay)/tth*0.5d0*dot_product(Nc, totstrain)
            
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

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_strainEnergy'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_strainEnergy
