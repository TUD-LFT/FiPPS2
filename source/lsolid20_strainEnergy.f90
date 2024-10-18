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
!> Calculates elementwise total strain energy for layered 20-node solid element Lsolid20
!
!> @details
!> Elementwise calculation of the total strain energy (tse) for layered 20-node solid
!> elements. The total strain energy is calculated as the strain tensor multiplied with
!> the stress tensor, multiplied by 0.5 and integrated over the element volume. Integration
!> over the element volume is performed numerically by using the Simpson rule out of plane
!> and Gauss rule in plane.
!
!> @author   Florian Dexl, TU Dresden, wiss. Mitarbeiter, 16.10.2023
!
!> $Id: lsolid20_stiff.f90 381 2018-11-01 15:05:25Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 381 $
!> $Date: 2018-11-01 16:05:25 +0100 (Do, 01. Nov 2018) $
!
! =================================================================================================
subroutine lsolid20_strainEnergy(fesim,nlay,Clay,lth,tth,node_coords,nip,nop,cid,displ,elem,scloop,lsolid20_tse)
! =================================================================================================
! use
!
 use integration_schemes
 use fesimulation_typen
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
! Data types
!
! Input
!
type(fe_simulation), intent(in)                   :: fesim
integer, intent(in)                               :: nlay !< Total number of layers in current element
integer, intent(in)                               :: nip !< Number of integration points in plane per direction
integer, intent(in)                               :: nop(nlay) !< (nlay)-Array with number of integration points over layer thickness for each layer
integer, intent(in)                               :: cid !< ID of coordinate system for orientation of element coosy
double precision, intent(in)                      :: Clay(6,6,nlay) !< (6,6,nlay)-Array which contains the material stiffness matrix in the element coordinate system for each layer
double precision, intent(in)                      :: lth(nlay) !< (nlay)-Array which contains the thickness for each layer
double precision, intent(in)                      :: tth !< Total thickness of laminat (sum of all values in lth)
double precision, intent(in)                      :: node_coords(20,3) !< (20,3)-Array containing global coordinates of element nodes
integer, intent(in)                               :: elem  !< Number of current element
double precision, intent(in), dimension(20*3)     :: displ !< Array containing the displacements in the nodal coordinate system
integer, intent(in)                               :: scloop
!
! Output
!
double precision, intent(out)                     :: lsolid20_tse !< Total strain energy of lsolid20 element
!
! Internal
!
integer                                           :: err_code=0, N, ii
integer                                           :: lay
double precision                                  :: fac
double precision, dimension(6,60)                 :: Bmat
double precision, allocatable, dimension(:)       :: xi_ip, eta_ip, zeta_ip
double precision, allocatable, dimension(:)       :: w_xi, w_eta, w_zeta
double precision, dimension(6)                    :: eps_gp_glob, stress_gp, eps_gp
double precision, dimension(6)                    :: initstress_elem, initstress_glob
double precision                                  :: weight,detJac
double precision, dimension(3,3)                  :: J
double precision                                  :: zeta_elem, lth_sum
!
! =================================================================================================
!
! Initialisation
!
 lsolid20_tse = 0.d0
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
     call lsolid20_strain_displ_matrix(xi_ip(ii),eta_ip(ii),zeta_elem,node_coords,Bmat,J(:,:),detjac)
            
     if (detJac < 0.d0) then
       err_code = 2
       goto 9999
     end if
     
     ! get strains at gauss-point in global coordinate system
     eps_gp_glob(:) = matmul(Bmat,displ)
     ! transform strains to element coordinate system
     call transform_eg_strain(fesim,eps_gp_glob(:),J(:,:),cid,eps_gp(:),'ge')

     ! Calculate stress vector in element coordinate system
     stress_gp(:) = matmul(Clay(:,:,lay),eps_gp(:))
      
     if ((fesim%is_multistep .EQ. .TRUE.) .AND. (fesim%lasten%subcases(scloop)%upstress .EQ. .TRUE.)) then
       ! If initial stress shall be applied
       ! Add initial stresses from last step
       ! Get initial stress in global coordinate system
       initstress_glob = fesim%residuals%initstr20s(elem)%lay(lay)%is_ip(:,ii)
       ! Transform initial stresses to element coordinate system
       call transform_eg_stress(fesim,initstress_glob,J(:,:),cid,initstress_elem,'ge')
       ! Add initial stress to stresses
       stress_gp(:) = stress_gp(:) + initstress_elem(:)
     end if

     fac = weight*detJac*lth(lay)/tth         
     
     lsolid20_tse = lsolid20_tse + 0.5d0*dot_product(stress_gp, eps_gp)*fac

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
   write(*,*)                      'lsolid20_strainEnergy'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!   
end if
!
return
!
end subroutine lsolid20_strainEnergy
