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
!> Computes strains and mechanical stress of 20-nodes layered solid element
!
!> @details
!> Computation of total stress (mechanical + thermal) strains, mechanical and
!> thermal strains and mechanical stress; the standard coordinate system for
!> all outputs is the element coordinate system, optionally the global coordinat
!> system is possible (chosen by logical globOut); strains and stresses of several
!> elements at a common node are averaged; all nodal values are calculated for the
!> element corner nodes; values at midside nodes are averaged from the values at
!> element corners; all values are calculated at 2x2 gauss points in the surface and
!> extrapolated to the element's corner nodes;
!> element values are the mean values of the corner node values;
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_strains_stresses.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_strains_stresses(fesim,dispg,nodal_temperatures,lsolid20temperatures,strain_me,strain_th,strain,stress,princStrain,princStress,strain_me_elem,strain_th_elem,strain_elem,stress_elem,princStrain_elem,princStress_elem,scloop)
! =================================================================================================
!
! use
 use konstanten
 use fesimulation_typen
 use integration_schemes
!
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Interface
!
  interface
    subroutine lsolid20_strain_displ_matrix(xi,eta,zeta,node_coords,Bmat,J,detjac)
      double precision, intent(in)                   :: xi, eta, zeta, node_coords(20,3)
      double precision, intent(out), dimension(6,60) :: Bmat
      double precision, intent(out), optional        :: detjac, J(3,3)
    end subroutine lsolid20_strain_displ_matrix
  end interface
!
! =================================================================================================
!
! Include
!
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation)                                            :: fesim
double precision, intent(in), dimension(fesim%num_dof)         :: dispg !< Array containing displacement of all nodes with respect to global coordinate system
double precision, dimension(fesim%num_nodes), intent(in)       :: nodal_temperatures !< Array containing nodal temperature of all nodes
double precision, dimension(size(fesim%elemente%lsolid20s,1)), intent(in)    :: lsolid20temperatures !< Array containing nodal temperature of all elements
integer, intent(in)                                            :: scloop
!
! Output
!
double precision, intent(out), dimension(fesim%num_nodes,6)    :: strain !< Total mechanical + thermal strain with respect to output coordinate system at all nodes; strain = B*u
double precision, intent(out), dimension(fesim%num_nodes,6)    :: strain_me !< Mechanical strain with respect to output coordinate system at all nodes; strain_me = strain - strain_th
double precision, intent(out), dimension(fesim%num_nodes,6)    :: strain_th !< Thermal strain with respect to output coordinate system at all nodes; strain_th = (th. exp. coeff)*(nodal temperature)
double precision, intent(out), dimension(fesim%num_nodes,6)    :: stress !< Mechanical stress with respect to output coordinate system at all nodes; stress = C*strain_me + initstress
double precision, intent(out), dimension(fesim%num_elements,6) :: strain_elem !< Element total strain with respect to output coordinate system
double precision, intent(out), dimension(fesim%num_elements,6) :: strain_me_elem !< Element mechanical strain with respect to output coordinate system
double precision, intent(out), dimension(fesim%num_elements,6) :: strain_th_elem !< Element thermal strains with respect to output coordinate system
double precision, intent(out), dimension(fesim%num_elements,6) :: stress_elem !< Element mechanical stress with respect to output coordinate system
double precision, intent(out), dimension(fesim%num_nodes,4)    :: princStrain !< Principal strains at all nodes
double precision, intent(out), dimension(fesim%num_nodes,4)    :: princStress !< Principal stresses at all nodes
double precision, intent(out), dimension(fesim%num_elements,4) :: princStrain_elem !< Principal strains at all elements
double precision, intent(out), dimension(fesim%num_elements,4) :: princStress_elem !< Principal stresses at all elements
!
! Internal
!
integer                                         :: err_code=0
integer                                         :: elem, node, ind, lay
integer                                         :: n_1, n_2, m_1, m_2, ii, jj
integer                                         :: j1, j2, nlay, node_id
integer                                         :: nip, nop, n_int, cid
integer                                         :: resLay, botLay, topLay
integer, dimension(fesim%num_nodes)             :: nodeNumbers
double precision                                :: detjac, zeta_elem, tth, lth_sum
double precision, dimension(6)                  :: eps_gp_glob, stress_gp_glob
double precision, dimension(6)                  :: initstress_elem, initstress_glob
double precision, dimension(20)                 :: ntemps
double precision, dimension(60)                 :: displ
double precision, allocatable, dimension(:,:,:) :: J
double precision, dimension(:,:), allocatable   :: eps_gp, eps_th_gp, eps_me_gp, stress_gp
double precision, dimension(6,20)               :: eps_node, eps_th_node, eps_me_node, stress_node
double precision, dimension(6,60)               :: Bmat
double precision, dimension(20,3)               :: node_coords
double precision, allocatable, dimension(:,:)   :: ath_lay
double precision, allocatable, dimension(:,:,:) :: Clay
double precision, allocatable, dimension(:)     :: lth, angles
double precision, allocatable, dimension(:)     :: xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta
logical                                         :: globOut
character(4)                                    :: outCS
!
! =================================================================================================
!
! Initialisation
!
 nip    = 2     ! number of integration points in plane
!
! =================================================================================================
!
! Calculation
 nodeNumbers = 0
 strain      = 0.d0
 strain_me   = 0.d0
 strain_th   = 0.d0
 stress      = 0.d0

 ! Loop over all Lsolid20 elements
 do elem = 1,size(fesim%elemente%lsolid20s,1)

   ! get number of layers for current element
   nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%lay

   if (nip .NE. 2) then
     write(*,*) 'wrong input on parameter nip (number of'
     write(*,*) 'integration points in plane)'
     write(*,*) 'stress calculation currently allows only'
     write(*,*) 'value 2 to be chosen'
     err_code = 2
     goto 9999
   end if
   
   ! allocate arrays dependent on layer number
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
   ! get output coordinate system
   globOut = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%globOut
   ! get number of result layer
   resLay  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%resLay
   
   ! get ID of bottom and top layer for output
   if (resLay .EQ. 0) then
     ! standard output displays top of top layer and bottom
     ! of bottom layer
     botLay = 1
     topLay = nlay
     
     ! output in element coordinate system
     outCS = 'elem'
     
   else if ((resLay .GE. 1) .AND. (resLay .LE. nlay)) then
     ! if specific layer number is given, top and bottom
     ! of this layer is chosen for output
     botLay = resLay
     topLay = resLay
     
     ! output in material coordinate system
     outCS = 'matr'

   else if (resLay .GT. nlay) then
     write(*,*) 'Invalid input for resLay'
     write(*,*) 'Layer number for output is greater than'
     write(*,*) 'number of defined layers'
     write(*,*) 'Please check input'
     err_code = 2
     goto 9999
   else if (resLay .LT. 0) then
     write(*,*) 'Invalid input for resLay'
     write(*,*) 'Layer number must be positive'
     write(*,*) 'Please check input'
     err_code = 2
     goto 9999
   end if
   
   ! check if output in global coordinate system is
   ! specified
   if (globOut .EQ. .TRUE.) then
     outCS = 'glob'
   end if

   ! Get Temperature at element nodes
   if (fesim%is_lsolid20temp == .true.) then
     ntemps(:)      = lsolid20temperatures(elem)
   else if (fesim%is_temperature == .true.) then
     ntemps = 0.d0
     do node = 1,20
       ntemps(node) = nodal_temperatures(fesim%elemente%lsolid20s(elem)%nids(node))
     end do
   end if
   
   ! get global displacement vector for element nodes
   do node = 1,20
     ind = (fesim%elemente%lsolid20s(elem)%nids(node)-1)*ndof
     do ii = 1,3
       displ(ii+(node-1)*3) = dispg(ind+ii)
     end do
     ! get node coordinates
     node_coords(node,:) = fesim%knoten%nodes(fesim%elemente%lsolid20s(elem)%nids(node))%coords(1:3)
   end do

   ! Loop over layers of current element
   lth_sum = 0.d0
   do lay = 1,nlay

     ! calculate summed thickness until current layer (needs to be processed for every layer!)
     lth_sum = lth_sum + lth(lay)
   
     ! if stress values inside the element are not needed, skip layers between top and bottom layer
     if ((fesim%is_multistep .EQ. .FALSE.) .AND. (.NOT. ((lay .EQ. botLay) .OR. (lay .EQ. topLay)))) then
       cycle
     end if

     ! get number of integration points over thickness for current layer
     nop  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(elem)%int_pid)%nop(lay)

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
     
     if (fesim%is_multistep .EQ. .TRUE.) then
       
       if (fesim%lasten%subcases(scloop)%upstress .EQ. .TRUE.) then ! if initial stress shall be applied
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
       
       ! Store stresses from current step at gauss points
       do ii = 1,n_int
         ! Transform stresses to global coordinate system
         call transform_eg_stress(fesim,stress_gp(:,ii),J(:,:,ii),cid,stress_gp_glob,'eg')
         ! Store with respect to global coordinate system
         fesim%residuals%initstr20s(elem)%lay(lay)%is_ip(:,ii) = stress_gp_glob(:)
       end do  
     end if

     ! Extrapolate strains and stresses from gauss points to corner nodes
     ! at the element top and bottom surface for output
     ! If specified, values are transformed to the global coordinate
     ! system
     if (lay .eq. botLay) then   ! bottom layer 
       n_1 = 1
       n_2 = nip*nip
       m_1 = 1
       m_2 = 4
       zeta_elem = -1.d0
       call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,eps_gp(:,n_1:n_2),   eps_node(:,m_1:m_2)   )
       call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,eps_th_gp(:,n_1:n_2),eps_th_node(:,m_1:m_2))
       call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,eps_me_gp(:,n_1:n_2),eps_me_node(:,m_1:m_2))
       call lsolid20_extrapolate_gp_stress(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,stress_gp(:,n_1:n_2),stress_node(:,m_1:m_2))
     end if
     if (lay .eq. topLay) then   ! top layer
       n_1 = nip*nip*(nop-1)+1
       n_2 = nip*nip*nop
       m_1 = 5
       m_2 = 8
       zeta_elem = 1.d0
       call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,eps_gp(:,n_1:n_2),   eps_node(:,m_1:m_2)   )
       call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,eps_th_gp(:,n_1:n_2),eps_th_node(:,m_1:m_2))
       call lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,eps_me_gp(:,n_1:n_2),eps_me_node(:,m_1:m_2))
       call lsolid20_extrapolate_gp_stress(fesim,node_coords,zeta_elem,cid,angles(lay),outCS,stress_gp(:,n_1:n_2),stress_node(:,m_1:m_2))
     end if
     
     ! calculate element values as mean of unaveraged nodal values
     do ii = 1,6
       strain_elem(fesim%elemente%lsolid20s(elem)%eid,ii)    = sum(   eps_node(ii,1:8))/8.d0
       strain_th_elem(fesim%elemente%lsolid20s(elem)%eid,ii) = sum(eps_th_node(ii,1:8))/8.d0
       strain_me_elem(fesim%elemente%lsolid20s(elem)%eid,ii) = sum(eps_me_node(ii,1:8))/8.d0
       stress_elem(fesim%elemente%lsolid20s(elem)%eid,ii)    = sum(stress_node(ii,1:8))/8.d0
     end do

     deallocate(eps_gp,eps_th_gp,eps_me_gp,stress_gp)
     deallocate(xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta)
     deallocate(J)

   ! End loop over layers
   end do
   
   ! Interpolate strains and stresses at mid nodes 9-20
   do jj=9,20
     if (jj .LE. 16) then
       ! Midside nodes in plane
       j1 = jj-8
       j2 = jj-7
       if (MOD((j2-1),4) == 0) j2 = j2-4
     else
       ! Midside nodes in thickness direction
       j1 = jj-16
       j2 = jj-12
     end if
       eps_node(:,jj)    = (eps_node(:,j1)    + eps_node(:,j2)   )/2.d0
       eps_th_node(:,jj) = (eps_th_node(:,j1) + eps_th_node(:,j2))/2.d0
       eps_me_node(:,jj) = (eps_me_node(:,j1) + eps_me_node(:,j2))/2.d0
       stress_node(:,jj) = (stress_node(:,j1) + stress_node(:,j2))/2.d0
   end do
   
   ! Sort element nodal strains and stresses
   do jj = 1,20
     node_id = fesim%elemente%lsolid20s(elem)%nids(jj)
     nodeNumbers(node_id)   = nodeNumbers(node_id) + 1
     strain(node_id,1:6)    = strain(node_id,1:6) + eps_node(:,jj)
     strain_th(node_id,1:6) = strain_th(node_id,1:6) + eps_th_node(:,jj)
     strain_me(node_id,1:6) = strain_me(node_id,1:6) + eps_me_node(:,jj)
     stress(node_id,1:6)    = stress(node_id,1:6) + stress_node(:,jj)
   end do
   
   ! deallocate dependent arrays
   deallocate(ath_lay,Clay,lth,angles)
 
 end do
 
 ! Average values at nodes
 do ii = 1,size(nodeNumbers,1)
   if (nodeNumbers(ii) == 0) nodeNumbers(ii) = 1
 end do
 !
 do ii = 1,6
   strain(:,ii)    = strain(:,ii)    / nodeNumbers
   strain_th(:,ii) = strain_th(:,ii) / nodeNumbers
   strain_me(:,ii) = strain_me(:,ii) / nodeNumbers
   stress(:,ii)    = stress(:,ii)    / nodeNumbers
 end do
 
 ! Calculate principal nodal strains and stresses
 do ii = 1, fesim%num_nodes
   call principalValues(strain_me(ii,:),stress(ii,:),princStrain(ii,:),princStress(ii,:))
 end do
 
 ! Calculate principal elemental strains and stresses
 do elem = 1, size(fesim%elemente%lsolid20s,1)
   call principalValues(strain_me_elem(fesim%elemente%lsolid20s(elem)%eid,:), &
                      & stress_elem(fesim%elemente%lsolid20s(elem)%eid,:), &
                      & princStrain_elem(fesim%elemente%lsolid20s(elem)%eid,:), &
                      & princStress_elem(fesim%elemente%lsolid20s(elem)%eid,:))
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_strains_stresses'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_strains_stresses
