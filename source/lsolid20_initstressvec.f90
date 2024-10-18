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
!> Calculates equivalent load vector for initial stresses
!
!> @details
!> Calculation of the equivalent load vector for consideration of given intial
!> stresses; initial stresses are taken from the global initstr20s, where the full
!> stress vector s = [s1,s2,s3,s23,s13,s12]^T is stored for every integration point
!> of every layer for every element; this array is filled during the previous
!> calculation step in subroutine lsolid20_strains_stresses; the integration of the
!> loadvector is done layerwise using the same integrationpoints as used for the
!> integration of the element stiffness matrix (see lsolid20_stiff);
!> zeta = [-1,1] is substituted for each layer such that the same integration points
!> can be used in every layer (see also lsolid20_stiff);
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_initstressvec.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_initstressvec(fesim,eid,nlay,lth,tth,node_coords,nip,nop,initstrvec)
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
INTERFACE 
  SUBROUTINE LSOLID20_STRAIN_DISPL_MATRIX(XI,ETA,ZETA,        &
               &NODE_COORDS,BMAT,J,DETJAC)
    REAL(KIND=8), INTENT(IN) :: XI
    REAL(KIND=8), INTENT(IN) :: ETA
    REAL(KIND=8), INTENT(IN) :: ZETA
    REAL(KIND=8), INTENT(IN) :: NODE_COORDS(20,3)
    REAL(KIND=8), INTENT(OUT) :: BMAT(6,60)
    REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: J(3,3)
    REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: DETJAC
  END SUBROUTINE LSOLID20_STRAIN_DISPL_MATRIX
END INTERFACE 
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
type(fe_simulation)                            :: fesim
integer, intent(in)                            :: eid !< Element ID of current element
integer, intent(in)                            :: nlay !< Total number of layers in current element
integer, intent(in)                            :: nip !< Number of integration points in plane per direction
integer, intent(in)                            :: nop(nlay) !< (nlay)-Array with number of integration points over layer thickness for each layer
double precision, intent(in)                   :: lth(nlay) !< (nlay)-Array which contains the thickness for each layer
double precision, intent(in)                   :: tth !< Total thickness of laminat (sum of all values in lth)
double precision, intent(in)                   :: node_coords(20,3) !< (20,3)-Array containing global coordinates of element nodes
! 
! Output
double precision, intent(out), dimension(60)   :: initstrvec !< Equivalent load vector for initial stresses of current element
!
! Internal
integer                                        :: lay, N
double precision, dimension(6)                 :: initstress
double precision, allocatable, dimension(:)    :: xi_ip, eta_ip, zeta_ip
double precision, allocatable, dimension(:)    :: w_xi ,w_eta, w_zeta
double precision, dimension(3,3)               :: J
double precision, dimension(6,60)              :: Bmat
double precision                               :: weight, detJac, lth_sum, zeta_elem

integer                                        :: err_code=0
integer                                        :: ii
!
! =================================================================================================
!
! Initialisation
!
 initstrvec = 0.d0
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

     ! get initial stresses for current integration point with
     ! respect to global coordinate system
     initstress(:) = fesim%residuals%initstr20s(eid)%lay(lay)%is_ip(:,ii)

     ! transform layerwise zeta-value of integration point zeta_ip(ii)
     ! to elementwise zeta_value zeta_elem
     ! (see eq.(22) in Panda, S.; Natarajan, R.: Analysis of laminated composite shell
     !  structures by finite element method. Computers & Structures 14(3-4). 1981. pp. 225-230)
     zeta_elem = -1.d0 + (2.d0*lth_sum - lth(lay)*(1.d0-zeta_ip(ii)))/tth

     ! Get strain-displacement matrix with respect to global coordinate system
     call lsolid20_strain_displ_matrix(xi_ip(ii),eta_ip(ii),zeta_elem,node_coords,Bmat,J,detJac)
     
     initstrvec = initstrvec + matmul(transpose(Bmat),initstress)*weight*detJac*lth(lay)/tth
     
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

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_initstressvec'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_initstressvec
