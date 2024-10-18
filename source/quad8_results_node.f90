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
!> Calculation of the nodal strains and stresses for the layered 8-Node Shell Element
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 29.06.2010
!
!> $Id: quad8_results_node.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_results_node (fesim, dispg, nodal_temperatures, quad8temperatures, strain, stress, pos, iLay)
! =================================================================================================
! use
!
use konstanten
use fesimulation_typen
use mat_func
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Input
!
  type(fe_simulation)                                        :: fesim
  double precision, dimension(fesim%num_dof), intent(in)     :: dispg                  !< (num_dof)-Array containing the nodal displacements in the nodal coordinate system
  double precision, dimension(fesim%num_nodes), intent(in)   :: nodal_temperatures     !< (num_nodes)-Array containing the nodal temperatures
  double precision, dimension(size(fesim%elemente%quad8s,1)), intent(in):: quad8temperatures !< (num_elements)-Array containing the elemental temperatures
  integer, intent(in)                                        :: pos                    !< Position flag for the results (1-TOP, 2-BOT, else-MID) 
  integer, intent(in)                                        :: iLay                   !< Layer number for which the results are calculated (iLay > 0 -> Layer results, iLay=0 -> if pos=1 then layer 1 BOT, if pos=2 then layer n TOP)
!
! Output
! 
  double precision, dimension(fesim%num_nodes,6), intent(out):: strain                !< (num_nodes,6)-Array containing the nodal strains (1- x-strain; 2- y-strain; 3- xy-strain; 4- 1st principal strain; 5- 2nd principal strain; 6- strain intensity)
  double precision, dimension(fesim%num_nodes,6), intent(out):: stress                !< (num_nodes,6)-Array containing the nodal stresses (1- x-stress; 2- y-stress; 3- xy-stress; 4- 1st principal stress; 5- 2nd principal stress; 6- stress intensity)
!
! Input + Output
!

!
! inner
!
  integer                                         :: elem, node
  integer                                         :: kk, ll, ind, numIP, layNum
  integer                                         :: err_code
  
  integer, dimension(:), allocatable              :: nodenumbers
  

  double precision, dimension(:), allocatable     :: disp_quad8
  double precision, dimension(8,6)                :: eps
  double precision, dimension(8,6)                :: sig
  double precision, dimension(4)                  :: princEpsilon, princSigma
  double precision                                :: zeta
  double precision, dimension(8,3)                :: node_coords
  double precision, dimension(6,6)                :: CTemp, TTemp
  double precision, dimension(8)                  :: ntemps
  double precision                                :: E, nue
  double precision, dimension(6)                  :: athTemp
  double precision, dimension(:,:),allocatable    :: Transform
  double precision, dimension(:,:), allocatable   :: prestrain, prestress
  double precision, dimension(6)                  :: epsvec_resort, sigvec_resort
  double precision                                :: tth, lth, lth_sum, zeta_elem
  double precision                                :: nu12, nu21, nu13, nu31, nu23, nu32
  double precision                                :: ym11, ym22, ym33
  double precision                                :: sm12, sm23, sm13
  double precision                                :: det, div
  double precision, dimension(2,2)                :: partCInv

  err_code = 0
  allocate(disp_quad8(8*ndof))
  allocate(nodenumbers(fesim%num_nodes))
  nodenumbers = 0

  ! Loop over all flat (reduced) quad8 elements
  do elem=1,size(fesim%elemente%quad8s,1)
    
    ! get global displacement vector for element nodes

    do node=1,8
       
      ind = (fesim%elemente%quad8s(elem)%nids(node)-1)*ndof
      
      do ll=1,ndof
         disp_quad8(ll+(node-1)*ndof) = dispg(ind+ll)
      end do

    end do
    
    do kk = 1,8
      node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(kk))%coords(1:3)
    end do
    
    ntemps = 0.d0
    if (fesim%is_quad8temp == .true.) then
      ntemps(:)      = quad8temperatures(elem)
    else if (fesim%is_temperature == .true.) then
      do node = 1,8
        ntemps(node) = nodal_temperatures(fesim%elemente%quad8s(elem)%nids(node))
      end do
    end if
    
    if (fesim%elemente%quad8s(elem)%propType == 1) then
    
      tth = fesim%eigenschaften%pshells(fesim%elemente%quad8s(elem)%int_pid)%mt
   
      E   = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(elem)%int_pid)%intMat1ID)%ym
      nue = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(elem)%int_pid)%intMat1ID)%nu
      
      ym11 = E
      ym22 = E
      ym33 = E
      
      nu12 = nue
      nu21 = nue
      nu31 = nue
      nu13 = nue
      nu32 = nue
      nu23 = nue
      
      sm12 = E/2.d0/(1.d0+nue)
      sm23 = E/2.d0/(1.d0+nue)
      sm13 = E/2.d0/(1.d0+nue)
      
      athTemp(1:3) = fesim%materialien%mat1s(fesim%eigenschaften%pshells(fesim%elemente%quad8s(elem)%int_pid)%intMat1ID)%ath
      athTemp(4:6) = 0.d0
      
      TTemp = 0.d0
      do kk = 1,6
        TTemp(kk,kk) = 1.d0
      end do
        
      zeta_elem = 0.d0 ! Mitte
      if (pos .eq. 1) then !Oberseite
        zeta_elem =  1.d0
      else if (pos .eq. 2) then !Unterseite
        zeta_elem = -1.d0
      end if
    
    else if (fesim%elemente%quad8s(elem)%propType == 2) then
  
      zeta = 0.d0 ! Mitte
      if (pos .eq. 1) then !Oberseite
        zeta =  1.d0
      else if (pos .eq. 2) then !Unterseite
        zeta = -1.d0
      end if
      
      if (iLay .eq. 0) then ! im Element
        if (pos .eq. 1) then !Oberseite
          layNum = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lay
        else if (pos .eq. 2) then !Unterseite
          layNum = 1
        else
          write(*,*) 'No results for layered quad8 element at middle surface!'
          err_code=1
          goto 9999
        end if
      else if (iLay .gt. 1) then
        if (iLay .gt. fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lay) then
          write(*,*) 'No results for layered quad8 element at layer: ', iLay
          err_code=2
          goto 9999
        end if
        layNum = iLay
      end if
    
      tth = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%tth
      lth = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lth(layNum)
      lth_sum = SUM(fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lth(1:layNum))
      zeta_elem = (2.d0*lth_sum - lth - tth)/tth + lth/tth*zeta
            
      CTemp = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%C(:,:,layNum)
      
      call M22INV(CTemp(1:2,1:2), partCInv, det)
      
      ym11 = 1.d0 / partCInv(1,1)
      ym22 = 1.d0 / partCInv(2,2)
      ym33 = ym11
      
      nu12 = CTemp(1,2)/CTemp(2,2)
      nu21 = CTemp(1,2)/CTemp(1,1)
      nu31 = nu21
      nu13 = nu31*ym11/ym33
      nu32 = nu21
      nu23 = nu32*ym22/ym33
      
      sm12 = CTemp(4,4)
      sm23 = CTemp(5,5)
      sm13 = CTemp(6,6)
    
      TTemp = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%T(:,:,layNum)
      
      athTemp(1:6) = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%ath(:,layNum)
      athTemp(3) = 0.d0
      athTemp(5:6) = 0.d0
      
    end if
      
    numIP = 2*2

    ! Material stiffness matrix according to:
    ! http://homes.civil.aau.dk/lda/Continuum/material.pdf
    ! Damkilde, L. Aalborg Universitet. 2008
    div = 1.d0 - nu23*nu32 - nu12*nu21 - nu13*nu31 - nu12*nu23*nu31 - nu21*nu32*nu13

    CTemp = 0.d0
    CTemp(1,1) = ym11*(1.d0 - nu23*nu32) / div        ! C11
    CTemp(1,2) = ym11*(nu23*nu31 + nu21) / div        ! C12
    CTemp(2,2) = ym22*(1.d0 - nu13*nu31) / div        ! C22
    CTemp(1,3) = ym11*(nu21*nu32 + nu31) / div        ! C13
    CTemp(2,3) = ym22*(nu12*nu31 + nu32) / div        ! C23
    CTemp(3,3) = ym33*(1.d0 - nu12*nu21) / div        ! C33
    CTemp(4,4) = sm12                                 ! C44
    CTemp(5,5) = sm23                                 ! C55
    CTemp(6,6) = sm13                                 ! C66
    CTemp(2,1) = CTemp(1,2)                            ! C21
    CTemp(3,1) = CTemp(1,3)                            ! C31
    CTemp(3,2) = CTemp(2,3)                            ! C32
    
    ! Transformieren der Steifigkeitmatrix in das Elementkoordinatensystem
    CTemp = matmul(transpose(TTemp),matmul(CTemp,TTemp))
    
    allocate(prestrain(numIP,6), prestress(numIP,6))
    
    call quad8_ip_strain (node_coords, &
                        & disp_quad8, &
                        & zeta_elem, &
                        & 2, &
                        & tth, &
                        & fesim%elemente%quad8s(elem)%area, &
                        & CTemp(3,1), &
                        & CTemp(3,2), &
                        & CTemp(3,3), &
                        & CTemp(3,4), &
                        & ntemps, &
                        & athTemp, &
                        & prestrain, &
                        & fesim%is_temperature, &
                        & fesim%is_quad8temp)
    
    if (fesim%elemente%quad8s(elem)%propType == 2) then
      prestrain(:,5) = 0.d0
      prestrain(:,6) = 0.d0
    end if
    
    allocate(Transform(8,numIP))
    call quad8_result_extrapolation_transform (2, Transform)
  
    eps(1:8,1:6) = matmul(Transform,prestrain(1:numIP,1:6))
    deallocate(Transform)
    
!     CTemp(:,3) = 0.d0
!     CTemp(3,:) = 0.d0
    
    do kk = 1,8
        sig(kk,1:6) = MATMUL(CTemp,eps(kk,1:6))
    end do
    
    do node = 1,8
      nodeNumbers(fesim%elemente%quad8s(elem)%nids(node)) = nodeNumbers(fesim%elemente%quad8s(elem)%nids(node)) + 1 
      strain(fesim%elemente%quad8s(elem)%nids(node),1:6)  = strain(fesim%elemente%quad8s(elem)%nids(node),1:6) + eps(node,1:6)
      stress(fesim%elemente%quad8s(elem)%nids(node),1:6)  = stress(fesim%elemente%quad8s(elem)%nids(node),1:6) + sig(node,1:6)
    end do
    
    deallocate(prestrain, prestress)

  end do !Elementschleife
  
  do ll = 1,size(nodeNumbers,1)
    if (nodeNumbers(ll) == 0) nodeNumbers(ll) = 1
  end do
  
  do ll = 1,size(strain,2)
    strain(:,ll) = strain(:,ll) / nodeNumbers
    stress(:,ll) = stress(:,ll) / nodeNumbers
  end do
  
!   write(*,'(A85)') '    NODE     EPTOX        EPTOY        EPTOZ        EPTOXY       EPTOYZ       EPTOXZ'
!   do ll = 1,size(strain,1)
!     write(*,'(I8,X,6(E13.5E3))') ll, strain(ll,1), strain(ll,2), strain(ll,3), strain(ll,4), strain(ll,5), strain(ll,6)
!   end do
!   
!   write(*,'(A85)') '    NODE     SX           SY           SZ           SXY          SYZ          SXZ'
!   do ll = 1,size(strain,1)
!     write(*,'(I8,X,6(E13.5E3))') ll, stress(ll,1), stress(ll,2), stress(ll,3), stress(ll,4), stress(ll,5), stress(ll,6)
!   end do
  
  do ll = 1,size(strain,1)
      
    epsvec_resort(1:3) = strain(ll,1:3)
    epsvec_resort(4)   = strain(ll,5)
    epsvec_resort(5)   = strain(ll,6)
    epsvec_resort(6)   = strain(ll,4)
    
    sigvec_resort(1:3) = stress(ll,1:3)
    sigvec_resort(4)   = stress(ll,5)
    sigvec_resort(5)   = stress(ll,6)
    sigvec_resort(6)   = stress(ll,4)
        
    call principalValues(epsvec_resort,sigvec_resort,princEpsilon,princSigma)
    
    strain(ll,1) = strain(ll,1)
    strain(ll,2) = strain(ll,2)
    strain(ll,3) = strain(ll,4)
    strain(ll,4) = princEpsilon(1)
    strain(ll,5) = princEpsilon(3)
    strain(ll,6) = princEpsilon(4)
    
    stress(ll,1) = stress(ll,1)
    stress(ll,2) = stress(ll,2)
    stress(ll,3) = stress(ll,4)
    stress(ll,4) = princSigma(1)
    stress(ll,5) = princSigma(3)
    stress(ll,6) = princSigma(4)
!    stress(ll,6) = MAX(ABS(princSigma(maxInd)-princSigma(midInd)), &
!                       ABS(princSigma(midInd)-princSigma(minInd)), &
!                       ABS(princSigma(minInd)-princSigma(maxInd)))
  end do
  
  deallocate(disp_quad8)
  deallocate(nodenumbers)
  
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_results_node'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return
  
end subroutine quad8_results_node
