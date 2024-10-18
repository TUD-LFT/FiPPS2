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
!> Calculation of the minimum reserve factor for the layered 8-Node Shell Element
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
!> $Id: quad8_results_elem_lay.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_results_elem_lay (fesim, dispg, nodal_temperatures, quad8temperatures, resFac, layNum, failTyp)
! =================================================================================================
! use
!
use konstanten
use fesimulation_typen
use mat_func
use failure_criteria
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Input
!
  type(fe_simulation)                                           :: fesim
  double precision, dimension(fesim%num_dof), intent(in)        :: dispg                  !< (num_dof)-Array containing the nodal displacements in the nodal coordinate system
  double precision, dimension(fesim%num_nodes), intent(in)      :: nodal_temperatures     !< (num_nodes)-Array containing the nodal temperatures
  double precision, dimension(size(fesim%elemente%quad8s,1)), intent(in)   :: quad8temperatures !< (num_nodes)-Array containing the nodal temperatures
!
! Output
! 
  double precision, dimension(fesim%num_elements), intent(out)  :: resFac                 !< (num_elements)-Array containing the minimum reserve factor for each element
  integer, dimension(fesim%num_elements), intent(out)           :: layNum                 !< (num_elements)-Array containing the layer number of the layer where the  minimum reserve factor was found for each element
  integer, dimension(fesim%num_elements), intent(out)           :: failTyp                !< (num_elements)-Array containing the failure where the minimum reserve factor was found for each element
!
! Input + Output
!

!
! inner
!
  integer                                         :: elem, node, numIP, iLay, pos
  integer                                         :: kk, ll, ind, typElemMin, layElemMin
  integer                                         :: typ
  integer                                         :: err_code
  

  double precision, dimension(:), allocatable     :: disp_quad8
  double precision, dimension(6)                  :: epsvec, epsvecLay, sigvec, sigvecLay
  double precision                                :: zeta
  double precision, dimension(8,3)                :: node_coords
  double precision, dimension(6,6)                :: CTemp, TTemp,CTT
  double precision, dimension(6)                  :: athTemp
  double precision, dimension(8)                  :: ntemps
  
  double precision, dimension(:,:), allocatable   :: prestrain
  double precision                                :: tth, lth, lth_sum, zeta_elem
  double precision                                :: nu12, nu21, nu13, nu31, nu23, nu32
  double precision                                :: ym11, ym22, ym33
  double precision                                :: sm12, sm23, sm13
  double precision                                :: det, div, rfElemMin
  double precision, dimension(2,2)                :: partCInv
  double precision, dimension(3)                  :: stresses, strains, globStrains
  double precision                                :: rf

  err_code = 0
  layNum = -1
  failTyp = -1
  resFac = -1.d0
  allocate(disp_quad8(8*ndof))

  ! Loop over all flat (reduced) quad8 elements
  do elem=1,size(fesim%elemente%quad8s,1)
  
    
    if (fesim%elemente%quad8s(elem)%propType == 1) then
      cycle
    end if
    
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
      ntemps(:) = quad8temperatures(elem)
    else if (fesim%is_temperature == .true.) then
      do node = 1,8
        ntemps(node) = nodal_temperatures(fesim%elemente%quad8s(elem)%nids(node))
      end do
    end if
    
    rfElemMin = 1.d300
    typElemMin = -1
    layElemMin = -1
    
    do iLay = 1, fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lay
        
      tth = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%tth
      lth = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lth(iLay)
      lth_sum = SUM(fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%lth(1:iLay))
            
      CTemp = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%C(:,:,iLay)
      
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
      
      athTemp(1:6) = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%ath(:,iLay)
      
      TTemp = fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%T(:,:,iLay)
        
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
      
      CTT = CTemp
!       CTT(:,3) = 0.d0
!       CTT(3,:) = 0.d0
    
      do pos = 1,2
      
        if (pos .eq. 1) then !Oberseite
          zeta =  1.d0
        else if (pos .eq. 2) then !Unterseite
          zeta = -1.d0
        end if
        
        zeta_elem = (2.d0*lth_sum - lth - tth)/tth + lth/tth*zeta
        
        allocate(prestrain(numIP,6))
        
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
        
        prestrain(:,5) = 0.d0
        prestrain(:,6) = 0.d0
        
     
        ! Laut 
        ! Hinton, E. ; Campbell, J. S.: Local and global smoothing of discontinuos finite
        ! element functions using a least squares method. In: international Journal for
        ! Numerical methods in Engineering 8 (1974), S. 461.
        ! ist zumindest für eine 2x2 Integration die Auswertung der geklätteten Funktion am
        ! Elementmittelpunkt gleich dem Mittelwert der Gauss-Punktwerte.

        epsvec(1) = SUM(prestrain(1:numIP,1))/numIP
        epsvec(2) = SUM(prestrain(1:numIP,2))/numIP
        epsvec(3) = SUM(prestrain(1:numIP,3))/numIP
        epsvec(4) = SUM(prestrain(1:numIP,4))/numIP
        epsvec(5) = SUM(prestrain(1:numIP,5))/numIP
        epsvec(6) = SUM(prestrain(1:numIP,6))/numIP
        
        globStrains(1) = epsvec(1)
        globStrains(2) = epsvec(2)
        globStrains(3) = epsvec(4)
        
        epsvecLay(1:6) = MATMUL(fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%T(1:6,1:6,iLay),epsvec(1:6))
        
        strains(1) = epsvecLay(1)
        strains(2) = epsvecLay(2)
        strains(3) = epsvecLay(4)
            
        ! Berechnung der Spannungen im Elementkoordinatensystem
        sigvec(1:6) = MATMUL(CTT,epsvec(1:6))
        
        sigvecLay(1:6) = MATMUL(transpose(fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%TInv(1:6,1:6,iLay)),sigVec(1:6))
        
        stresses(1) = sigvecLay(1)
        stresses(2) = sigvecLay(2)
        stresses(3) = sigvecLay(4)
        
        call getRF_mat8(stresses, &
                      & strains, &
                      & globStrains, &
                      & fesim%eigenschaften%pcomps(fesim%elemente%quad8s(elem)%int_pid)%intMatID(iLay), &
                      & rf, &
                      & typ, &
                      & fesim%versagenskriterien, &
                      & fesim%materialien%mat8s)
        
        if (rf .lt. rfElemMin) then
          rfElemMin = rf
          typElemMin = typ
          layElemMin = iLay
        end if
        
        deallocate(prestrain)
      
      end do !Positionsschleife
      
    end do !Lagenschleife
    
    resFac(fesim%elemente%quad8s(elem)%eid) = rfElemMin
    layNum(fesim%elemente%quad8s(elem)%eid) = layElemMin
    failTyp(fesim%elemente%quad8s(elem)%eid) = typElemMin

  end do !Elementschleife
  
  deallocate(disp_quad8)

  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_results_elem_lay'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return
  
end subroutine quad8_results_elem_lay
