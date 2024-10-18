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
!> Calculates element geometric stiffness matrix for layered 20-node solid element Lsolid20
!
!> @details
!> Calculation of element geometric stiffness matrix for layered 20-node solid element
!> Lsolid20; integration is done layerwise, zeta = [-1,1] is substituted for
!> each layer to run over thickness of layer such that the same integration point
!> values can be used in every layer; by this, the layer thicknesses are taken into
!> account as relative to the total laminat thickness which therefore needs not to be
!> equal to the actual thickness of the element as layer thicknesses are scaled to it;
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 27.09.2018
!
!> $Id: lsolid20_geostiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_geostiff(fesim,scloop,nlay,Clay,ath_lay,lth,tth,node_coords,nip,nop,cid,lsolid20_Kg,disp_l20,nodetemps,initstress_glob)
! =================================================================================================
! use
!
 use integration_schemes
 use fesimulation_typen
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
    subroutine lsolid20_strain_displ_matrix(xi,eta,zeta,node_coords,Gmat,J,detjac)
      double precision, intent(in)                   :: xi, eta, zeta, node_coords(20,3)
      double precision, intent(out), dimension(6,60) :: Gmat
      double precision, intent(out), optional        :: detjac, J(3,3)
    end subroutine lsolid20_strain_displ_matrix
  end interface
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation), intent(in)                   :: fesim             !< FE-Simulation data
integer, intent(in)                               :: scloop            !< Number of current subcase
integer, intent(in)                               :: nlay              !< Total number of layers in current element
integer, intent(in)                               :: nip               !< Number of integration points in plane per direction
integer, intent(in)                               :: nop(nlay)         !< (nlay)-Array with number of integration points over layer thickness for each layer
integer, intent(in)                               :: cid               !< ID of coordinate system for orientation of element coosy
double precision, intent(in)                      :: Clay(6,6,nlay)    !< (6,6,nlay)-Array which contains the material stiffness matrix in the element coordinate system for each layer
double precision, intent(in)                      :: ath_lay(6,nlay)   !< (6,nlay)-Array which contains the thermal expansion coefficients for each layer
double precision, intent(in)                      :: lth(nlay)         !< (nlay)-Array which contains the thickness for each layer
double precision, intent(in)                      :: tth               !< Total thickness of laminat (sum of all values in lth)
double precision, intent(in)                      :: node_coords(20,3) !< (20,3)-Array containing global coordinates of element nodes
double precision, dimension(60), intent(in)       :: disp_l20          !< global displacements for each node
double precision, dimension(20), intent(in)       :: nodetemps         !< (20)-Array containing the temperature of element nodes
double precision, dimension(nlay,6,nip*nip*maxval(nop)), intent(in) :: initstress_glob   !< Array containing the initial stresses for every integration point in each layer
!
! Output
!
double precision, intent(out), dimension(60,60)   :: lsolid20_Kg       !< 60x60-Element geometric stiffness matrix (no values for rotational dofs)
!
! Internal
!
integer                                           :: err_code=0, N, ii
integer                                           :: jj, kk
integer                                           :: lay
double precision, dimension(6)                    :: eps_gp, eps_th_gp, eps_me_gp, eps_gp_glob, stress_gp, stress_gp_glob, initstress_elem
double precision, allocatable, dimension(:)       :: xi_ip, eta_ip, zeta_ip
double precision, allocatable, dimension(:)       :: w_xi, w_eta, w_zeta
double precision, dimension(20)                   :: dNidx, dNidy, dNidz
double precision, dimension(3,3)                  :: J, lam
double precision, dimension(9,9)                  :: S
double precision, dimension(9,60)                 :: Gmat
double precision, dimension(6,60)                 :: StrainDisplMat
double precision                                  :: weight,detJac
double precision                                  :: zeta_elem, lth_sum
!
! =================================================================================================
!
! Initialisation
!
 lsolid20_Kg = 0.d0
 lth_sum = 0.d0
!
! =================================================================================================
!
! Calculation
!
 do lay = 1,nlay  

   ! Total number of integration points in current layer
   N = nip*nip*nop(lay)

   ! allocate layer specific arrays for integration point values
   allocate(xi_ip(N),eta_ip(N),zeta_ip(N),w_xi(N),w_eta(N),w_zeta(N))
  
   ! get integration points and corresponding weights
   call integration_points_3d(     nip, "Gauss",   xi_ip, w_xi, &
                                   nip, "Gauss",  eta_ip, w_eta, &
                              nop(lay), "Simps", zeta_ip, w_zeta)

   ! calculate summed thickness until current layer
   lth_sum = lth_sum + lth(lay)

   if (nip .NE. 2) then
     write(*,*) 'wrong input on parameter nip (number of'
     write(*,*) 'integration points in plane)'
     write(*,*) 'stress calculation currently allows only'
     write(*,*) 'value 2 to be chosen'
     err_code = 2
     goto 9999
   end if

   ! Loop over integration points
   do ii = 1,N

     weight = w_xi(ii)*w_eta(ii)*w_zeta(ii)

     ! transform layerwise zeta-value of integration point zeta_ip(ii)
     ! to elementwise zeta_value zeta_elem
     ! (see eq.(22) in Panda, S.; Natarajan, R.: Analysis of laminated composite shell
     !  structures by finite element method. Computers & Structures 14(3-4). 1981. pp. 225-230)
     zeta_elem = -1.d0 + (2.d0*lth_sum - lth(lay)*(1.d0-zeta_ip(ii)))/tth

     ! get strain-displacement matrix with respect to global coordinate system and jacobi matrix
     ! at current gauss point
     call lsolid20_strain_displ_matrix(xi_ip(ii),eta_ip(ii),zeta_elem,node_coords,StrainDisplMat,J(:,:),detjac)
     ! get strains at gauss-point in global coordinate system
     eps_gp_glob(:) = matmul(StrainDisplMat,disp_l20)
     ! transform strains to element coordinate system
     call transform_eg_strain(fesim,eps_gp_glob(:),J(:,:),cid,eps_gp(:),'ge')
     
     ! Calculation of thermal strains at gauss points of current surface
     if ((fesim%is_temperature == .true.) .or. (fesim%is_lsolid20temp == .true.)) then
      ! Calculate thermal strain in element coordinate system  
       call lsolid20_eps_th(nodetemps,xi_ip(ii),eta_ip(ii),zeta_elem,ath_lay(:,lay),eps_th_gp(:))
     else
       eps_th_gp(:) = 0.d0
     end if

     ! Stress calculation
     ! Calculate mechanical strains in element coordinate system
     eps_me_gp(:) = eps_gp(:) - eps_th_gp(:)
     ! Calculate stress vector in element coordinate system
     stress_gp(:) = matmul(Clay(:,:,lay),eps_me_gp(:))

     if (fesim%is_multistep .EQ. .TRUE.) then
       if (fesim%lasten%subcases(scloop)%upstress .EQ. .TRUE.) then ! if initial stress shall be applied
         ! Add initial stresses from last step
         ! Transform initial stresses to element coordinate system
         call transform_eg_stress(fesim,initstress_glob(lay,:,ii),J(:,:),cid,initstress_elem,'ge')
         ! Add initial stress to stresses
         stress_gp(:) = stress_gp(:) + initstress_elem(:)
       end if
     end if

     ! transform stresses to global coordinate system
     call transform_eg_stress(fesim,stress_gp(:),J(:,:),cid,stress_gp_glob(:),'eg')
     
     ! set up stress tensor
     lam(1,1) = stress_gp_glob(1); lam(1,2) = stress_gp_glob(6); lam(1,3) = stress_gp_glob(5)
     lam(2,1) = stress_gp_glob(6); lam(2,2) = stress_gp_glob(2); lam(2,3) = stress_gp_glob(4)
     lam(3,1) = stress_gp_glob(5); lam(3,2) = stress_gp_glob(4); lam(3,3) = stress_gp_glob(3)

     ! set up stress matrix
     S(:,:)     = 0.d0
     S(1:3,1:3) = lam(:,:)
     S(4:6,4:6) = lam(:,:)
     S(7:9,7:9) = lam(:,:)
  
     ! Get derivates of shape functions with respect to global coordinates
     call hexa20_shapefunctions_x_y_z(xi_ip(ii),eta_ip(ii),zeta_elem,node_coords,dNidx,dNidy,dNidz,detJac,J)

     ! set up G-Matrix
     Gmat(:,:) = 0.d0
     do kk = 1,3
       do jj = 1,20
         Gmat((kk-1)*3+1,(jj-1)*3+kk) = dNidx(jj)
         Gmat((kk-1)*3+2,(jj-1)*3+kk) = dNidy(jj)
         Gmat((kk-1)*3+3,(jj-1)*3+kk) = dNidz(jj)
       end do
     end do
     
     ! integrate geometric stiffness matrix
     lsolid20_Kg = lsolid20_Kg + matmul(matmul(transpose(Gmat),S),Gmat)*weight*detJac*lth(lay)/tth

   end do
   
   ! deallocate layer specific arrays for integration point values
   deallocate(xi_ip,eta_ip,zeta_ip,w_xi,w_eta,w_zeta)
   
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue
!
if (err_code /= 0) then
!   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_geostiff'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!   
end if
!
return
!
end subroutine lsolid20_geostiff
