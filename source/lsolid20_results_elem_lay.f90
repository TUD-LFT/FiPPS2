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
!> Calculation of the minimum reserve factor for the layered 20-Node Solid Element
!
!> @details
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 07.05.2017
!
!> $Id: lsolid20_results_elem_lay.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine  lsolid20_results_elem_lay (fesim, dispg, nodal_temperatures, lsolid20temperatures, resFac, layNum, failTyp, scloop)
! =================================================================================================
! use
!
use konstanten
use fesimulation_typen
use integration_schemes
use failure_criteria
!
! =================================================================================================
!
  implicit none
  
  interface 
    subroutine lsolid20_strain_displ_matrix(xi,eta,zeta,node_coords,Bmat,J,detjac)
      REAL(KIND=8), INTENT(IN) :: XI
      REAL(KIND=8), INTENT(IN) :: ETA
      REAL(KIND=8), INTENT(IN) :: ZETA
      REAL(KIND=8), INTENT(IN) :: NODE_COORDS(20,3)
      REAL(KIND=8), INTENT(OUT) :: BMAT(6,60)
      REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: J(3,3)
      REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: DETJAC
    end subroutine lsolid20_strain_displ_matrix
  end interface
!
! =================================================================================================
!
! Input
!
  type(fe_simulation)                                           :: fesim
  double precision, dimension(fesim%num_dof), intent(in)        :: dispg                  !< (num_dof)-Array containing the nodal displacements in the nodal coordinate system
  double precision, dimension(fesim%num_nodes), intent(in)      :: nodal_temperatures     !< (num_nodes)-Array containing the nodal temperatures
  double precision, dimension(fesim%num_elements), intent(in)   :: lsolid20temperatures   !< (num_elements)-Array containing the nodal temperatures
  integer, intent(in)                                           :: scloop
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
  integer                                         :: elem, node, numIP, lay, pos
  integer                                         :: ii, kk, ll, ind, typElemMin, layElemMin
  integer                                         :: nop, nip, n_int, n_1, n_2
  integer                                         :: nlay, cid, side
  integer                                         :: typ
  integer                                         :: err_code
  

  double precision, dimension(20*3)               :: displ
  double precision, dimension(6)                  :: epsvec, epsvecLay, sigvec, sigvecLay
  double precision, dimension(20,3)               :: node_coords
  double precision, dimension(20)                 :: ntemps
  
  double precision, dimension(:,:), allocatable   :: eps_gp,eps_th_gp,eps_me_gp,stress_gp
  double precision, dimension(:), allocatable     :: xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta
  double precision, allocatable, dimension(:,:,:) :: J
  double precision, allocatable, dimension(:,:)   :: ath_lay
  double precision, allocatable, dimension(:,:,:) :: Clay
  double precision, allocatable, dimension(:)     :: lth, angles
  double precision, dimension(6)                  :: eps_gp_glob, stress_gp_glob
  double precision, dimension(6)                  :: initstress_elem, initstress_glob
  double precision, dimension(6,4)                :: eps_me_node_mat, eps_me_node_elem, stress_node
  double precision, dimension(6,60)               :: Bmat
  double precision                                :: tth, lth_sum, zeta_elem
  double precision                                :: nu12, nu21, nu13, nu31, nu23, detJac, nu32
  double precision                                :: ym11, ym22, ym33
  double precision                                :: sm12, sm23, sm13
  double precision                                :: rfElemMin
  double precision, dimension(6)                  :: stresses, strains, globStrains
  double precision                                :: rf

  err_code = 0
  layNum = -1
  failTyp = -1
  resFac = -1.d0

  ! Loop over lsolid20 num_elements
  do elem=1,size(fesim%elemente%lsolid20s,1)

    ! get number of layers for current element
    nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%lay

    allocate(ath_lay(6,nlay))
    allocate(Clay(6,6,nlay))
    allocate(lth(nlay))
    allocate(angles(nlay))
    
    ! get thermal expansion coefficients for layers
    if ((fesim%is_temperature == .true.) .or. (fesim%is_lsolid20temp == .true.)) ath_lay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%ath(:,:)
    ! get stiffness matrices for layers
    Clay    = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%C(:,:,:)
    ! get thickness for layers
    lth     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%lth(:)
    ! get total thickness of laminat
    tth     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%tth
    ! get angles of layers
    angles  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%angle(:)
    ! get ID of coordinate system for orientation of element coosy
    cid     = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%cid

    ! get global displacement vector for element nodes
    do node=1,20
      ind = (fesim%elemente%lsolid20s(elem)%nids(node)-1)*ndof
      do ll=1,3
        displ(ll+(node-1)*3) = dispg(ind+ll)
      end do
      ! get node coordinates
      node_coords(node,:) = fesim%knoten%nodes(fesim%elemente%lsolid20s(elem)%nids(node))%coords(1:3)
    end do

    ! Get Temperature at element nodes
    ntemps = 0.d0
    if (fesim%is_lsolid20temp == .true.) then
      ntemps(:)      = lsolid20temperatures(elem)
    else if (fesim%is_temperature == .true.) then
      do node = 1,20
        ntemps(node) = nodal_temperatures(fesim%elemente%lsolid20s(elem)%nids(node))
      end do
    end if
    
    rfElemMin = 1.d300
    typElemMin = -1
    layElemMin = -1
    
    do lay = 1, fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%lay
    
      ! calculate summed thickness until current layer
      lth_sum = SUM(fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%lth(1:lay))

      ! get number of integration points over thickness for current layer
      nop  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%nop(lay)
      
      ! number of integration points in plane
      nip = 2

      n_int = nip*nip*nop
     
      allocate(eps_gp(6,n_int),eps_th_gp(6,n_int),eps_me_gp(6,n_int),stress_gp(6,n_int))
      allocate(xi_ip(n_int),eta_ip(n_int),zeta_ip(n_int),w_xi(n_int),w_eta(n_int),w_zeta(n_int))
      allocate(J(3,3,n_int))
      
      call integration_points_3d(nip, "Gauss",   xi_ip, w_xi, &
                                 nip, "Gauss",  eta_ip, w_eta, &
                                 nop, "Simps", zeta_ip, w_zeta)
      
      ! Loop over gauss points of current layer
      do ii = 1,n_int
        zeta_elem = -1.d0 + (2.d0*lth_sum - lth(lay)*(1.d0-zeta_ip(ii)))/tth
        ! get strain-displacement matrix with respect to global coordinate system and jacobi matrix
        ! at current gauss point
        call lsolid20_strain_displ_matrix(xi_ip(ii),eta_ip(ii),zeta_elem,node_coords,Bmat,J(:,:,ii),detjac)
        ! get strains at gauss-point in global coordinate system
        eps_gp_glob(:) = matmul(Bmat,displ)
        ! transform strains to element coordinate system
        call transform_eg_strain(fesim,eps_gp_glob(:),J(:,:,ii),cid,eps_gp(:,ii),'ge')
      end do
      
      ! Calculation of thermal strains at gauss points of current surface
      ! Loop over gauss points of current layer
      if ((fesim%is_temperature == .true.) .or. (fesim%is_lsolid20temp == .true.)) then  
        do ii = 1,n_int
          zeta_elem = -1.d0 + (2.d0*lth_sum - lth(lay)*(1.d0-zeta_ip(ii)))/tth
          ! Calculate thermal strain in element coordinate system  
          call lsolid20_eps_th(ntemps,xi_ip(ii),eta_ip(ii),zeta_elem,ath_lay(:,lay),eps_th_gp(:,ii))
        end do
      else
        eps_th_gp(:,:) = 0.d0
      end if
      
      ! Stress calculation
      ! Loop over gauss points of current layer
      do ii = 1,n_int
        ! Calculate mechanical strains in element coordinate system
        eps_me_gp(:,ii) = eps_gp(:,ii) - eps_th_gp(:,ii)
        ! Calculate stress vector in element coordinate system
        stress_gp(:,ii) = matmul(Clay(:,:,lay),eps_me_gp(:,ii))
      end do
      
      if ((fesim%is_multistep .EQ. .TRUE.) .AND. (fesim%lasten%subcases(scloop)%upstress .EQ. .TRUE.)) then
        ! If initial stress shall be applied
        ! Add initial stresses from last step
        do ii = 1,n_int
          ! Get initial stress in global coordinate system
          initstress_glob = fesim%residuals%initstr20s(elem)%lay(lay)%is_ip(:,ii)
          ! Transform initial stresses to element coordinate system
          call transform_eg_stress(fesim,initstress_glob,J(:,:,ii),cid,initstress_elem,'ge')
          ! Add initial stress to stresses
          stress_gp(:,ii) = stress_gp(:,ii) + initstress_elem(:)
        end do 
      end if
      
      ! Extrapolate strains and stresses from gauss points to corner nodes
      ! at the element top and bottom surface
      ! system
      
      do side=1,2
        if (side .eq. 1) then ! bottom layer
          n_1 = 1
          n_2 = nip*nip
          zeta_elem = -1.d0
        else if (side .eq. 2) then ! top layer
          n_1 = nip*nip*(nop-1)+1
          n_2 = nip*nip*nop
          zeta_elem = 1.d0
        end if

        call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),'matr',eps_me_gp(:,n_1:n_2),eps_me_node_mat(:,:))
        call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),'elem',eps_me_gp(:,n_1:n_2),eps_me_node_elem(:,:))
        call lsolid20_extrapolate_gp_stress(fesim,node_coords,zeta_elem,cid,angles(lay),'matr',stress_gp(:,n_1:n_2),stress_node(:,:))
      
        ! calculate element values as mean of unaveraged nodal values
        do ii = 1,6
          strains(ii)     = sum(eps_me_node_mat(ii,:))/4.d0
          globStrains(ii) = sum(eps_me_node_elem(ii,:))/4.d0
          stresses(ii)    = sum(stress_node(ii,:))/4.d0
        end do
        
        call getRF_mat20(stresses, &
                       & strains, &
                       & globStrains, &
                       & fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%intMatID(lay), &
                       & rf, &
                       & typ, &
                       & fesim%versagenskriterien, &
                       & fesim%materialien%mat20s)

        if (rf .lt. rfElemMin) then
          rfElemMin = rf
          typElemMin = typ
          layElemMin = lay
        end if
      end do
      
      deallocate(eps_gp,eps_th_gp,eps_me_gp,stress_gp)
      deallocate(xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta)
      deallocate(J)

    end do
    
    resFac(fesim%elemente%lsolid20s(elem)%eid)  = rfElemMin
    layNum(fesim%elemente%lsolid20s(elem)%eid)  = layElemMin
    failTyp(fesim%elemente%lsolid20s(elem)%eid) = typElemMin

    deallocate(ath_lay,Clay,lth,angles)
  end do

  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'lsolid20_results_elem_lay'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return
  
end subroutine lsolid20_results_elem_lay
