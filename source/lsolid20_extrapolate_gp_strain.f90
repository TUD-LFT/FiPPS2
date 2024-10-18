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
!> Extrapolate strains from 2x2 gauss points in a xi-eta-plane to the corner nodes
!> of the same plane; optional transformation to another coordinate system
!
!> @details
!> Strains are extrapolated from the 2x2 gauss points in a xi-eta-plane to the corner
!> nodes of the same plane; optionally, the resulting strains are transformed to the
!> element or material coordinate system; input strains must be given in the element
!> coordinate system and respect the following ordering:
!> s = [s11, s22, s33, s23, s13, s12]^T
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: lsolid20_extrapolate_gp_strain.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_extrapolate_gp_strain(fesim,node_coords,zeta,cid,phi,outCS,strains_gp,strains_cn)
! =================================================================================================
! use
!
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
! Include
!
!
! =================================================================================================
!
! Data types
!
! Input
!
type(fe_simulation), intent(in)                 :: fesim
double precision, intent(in), dimension(20,3)   :: node_coords !< (20,3)-array containing the coordinates of the element's nodes
double precision, intent(in)                    :: zeta !< zeta-coordinate of the xi-eta-plane (zeta = const.)
integer, intent(in)                             :: cid !< ID of coordinate system for orientation of element coosy
double precision                                :: phi !< Angle of layer (between material and coordinate system)
character(4)                                    :: outCS !< Specifier for output coordinate system: 'elem': element coordinate system, 'matr': layer (material) coordinate system, 'glob': global coordinate system
double precision, intent(in), dimension(6,4)    :: strains_gp !< (6,4)-array containing the strains at the plane's 2x2 gauss points with respect to the element coordinate system and the following ordering: s = [s11, s22, s33, s23, s13, s12]^T
!
! Output
!
double precision, intent(out), dimension(6,4)   :: strains_cn !< (6,4)-array containing the strains at the plane's corner nodes with respect to the specified output coordinate system
!
! Internal
!
integer                                         :: err_code=0
integer                                         :: ii
double precision                                :: detjac
double precision, dimension(20)                 :: dNidx, dNidy, dNidz
double precision, dimension(3,3)                :: E
double precision, dimension(4,4)                :: EXT
double precision, dimension(6,4)                :: strains_temp
double precision, dimension(6,6)                :: T
! double precision, parameter                     :: fac=1.d0/sqrt(3.d0)
! double precision, parameter, dimension(4)       :: xi_gp =(/ -fac, fac, fac,-fac /)
! double precision, parameter, dimension(4)       :: eta_gp=(/ -fac,-fac, fac, fac /)
double precision, dimension(4)                  :: xi_gp
double precision, dimension(4)                  :: eta_gp
double precision, dimension(4)                  :: wDummy
!
! =================================================================================================
!
! Initialisation
!
 ! Matrix for extrapolation from integration points to element nodes (2x2-Integration)
!  EXT(1,1) = 1.d0+sqrt(3.d0)/2.d0; EXT(1,2) = -0.5d0;               EXT(1,3) = 1.d0-sqrt(3.d0)/2.d0; EXT(1,4) = -0.5d0
!  EXT(2,1) = -0.5d0;               EXT(2,2) = 1.d0+sqrt(3.d0)/2.d0; EXT(2,3) = -0.5d0;               EXT(2,4) = 1.d0-sqrt(3.d0)/2.d0
!  EXT(3,1) = 1.d0-sqrt(3.d0)/2.d0; EXT(3,2) = -0.5d0;               EXT(3,3) = 1.d0+sqrt(3.d0)/2.d0; EXT(3,4) = -0.5d0
!  EXT(4,1) = -0.5d0;               EXT(4,2) = 1.d0-sqrt(3.d0)/2.d0; EXT(4,3) = -0.5d0;               EXT(4,4) = 1.d0+sqrt(3.d0)/2.d0
 EXT(1,1) = 1.d0+sqrt(3.d0)/2.d0; EXT(1,2) = -0.5d0;               EXT(1,3) = -0.5d0;               EXT(1,4) = 1.d0-sqrt(3.d0)/2.d0
 EXT(2,1) = -0.5d0;               EXT(2,2) = 1.d0+sqrt(3.d0)/2.d0; EXT(2,3) = 1.d0-sqrt(3.d0)/2.d0; EXT(2,4) = -0.5d0
 EXT(3,1) = 1.d0-sqrt(3.d0)/2.d0; EXT(3,2) = -0.5d0;               EXT(3,3) = -0.5d0;               EXT(3,4) = 1.d0+sqrt(3.d0)/2.d0
 EXT(4,1) = -0.5d0;               EXT(4,2) = 1.d0-sqrt(3.d0)/2.d0; EXT(4,3) = 1.d0+sqrt(3.d0)/2.d0; EXT(4,4) = -0.5d0
 
 call integration_points_2d(2, "Gauss",  xi_gp, wDummy, &
                            2, "Gauss", eta_gp, wDummy)
!
! =================================================================================================
!
! Calculation

 ! extrapolate strains from gauss-points to element corner nodes
 do ii = 1,6
   strains_temp(ii,:) = matmul(EXT,strains_gp(ii,:))
 end do

 if (outCS .EQ. 'elem') then
 
   ! strains are already in element coordinate system
   strains_cn(:,:) = strains_temp(:,:)
   
 else if (outCS .EQ. 'matr') then
 
   ! get cosine matrix
   E      = 0.d0
   E(1,1) =  cos(phi); E(1,2) = sin(phi);
   E(2,1) = -sin(phi); E(2,2) = cos(phi);
   E(3,3) = 1.d0
   
   ! get transformation matrix
   call lsolid20_transformation_matrix(E,T)
   
   ! transform strain
   strains_cn = matmul(T,strains_temp)
   
 else if (outCS .EQ. 'glob') then
   ! ---------------------------------------------------------------------------
   ! Begin Retransformation of strains from element to global coordinate system
   ! loop over gauss points in current surface
   do ii = 1,4
     ! Get jacobi matrix at gauss point
     call hexa20_shapefunctions_x_y_z(xi_gp(ii),eta_gp(ii),zeta,node_coords,dNidx,dNidy,dNidz,detjac,E)
     ! Transform vector of strains from
     ! element to global coordinate system
     call transform_eg_strain(fesim,strains_temp(:,ii),E,cid,strains_cn(:,ii),'eg')
   end do
   ! End Transformation of strains from element to global coordinate system
   ! ---------------------------------------------------------------------------
 else
   write(*,*) 'wrong specifier for output coordinate system'
   write(*,*) 'valid specifiers are:'
   write(*,*) 'outCS = ''elem'': element coordinate system'
   write(*,*) 'outCS = ''matr'': layer (material) coordinate system'
   write(*,*) 'outCS = ''glob'': global coordinate system'
   err_code = 2
   goto 9999
 end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_extrapolate_gp_strain'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_extrapolate_gp_strain
