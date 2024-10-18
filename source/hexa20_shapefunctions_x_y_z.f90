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
!> Subroutine for computing the shapefunctions of a 20-node hexaedron.
!
!> @details
!> Computes the derivations of the shapefunctions with respect to the global
!> coordinates. By the way, the jacobian determinant is computed for 
!> checking the element geometry.
!
!> @author Thomas Dylla, TU Dresden, Diplomarbeit, 13.12.2013
!
!> @author Florian Dexl, (added output of Jacobi-matrix), 2015
!
!> $Id: hexa20_shapefunctions_x_y_z.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine hexa20_shapefunctions_x_y_z(xi,eta,zeta,node_coordinates,dNidx,dNidy,dNidz,detjac,jac)
! =================================================================================================
! use
!
use globale_variablen
use mat_func
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Interface
!   This interface has to be implemented for calling solid20_shapefunctions.
!   This subroutine has optional arguments. It has to be defined for checking implemented data types.    
interface
    subroutine hexa20_shapefunctions_xi_eta_zeta(xi, eta, zeta, Ni, dNidxi, dNideta, dNidzeta)
      double precision, intent(in)                  :: xi, eta, zeta
      double precision, dimension(20), optional, intent(out)        :: Ni
      double precision, dimension(20), optional, intent(out)        :: dNidxi, dNideta, dNidzeta
    end subroutine hexa20_shapefunctions_xi_eta_zeta
end interface   
!
! =================================================================================================
!
! Data types
!
! Input
!
double precision                      :: xi !< Natural coordinate xi
double precision                      :: eta !< Natural coordinate eta
double precision                      :: zeta !< Natural coordinate zeta
double precision,dimension(20,3)              :: node_coordinates !<  Global coordinates of the nodes assigned to the element. Col 1: x ; Col 2: y ; Col 3: z
!
! Output
!
double precision, intent(out),dimension(20)       :: dNidx !< Derivations of the shapefunctions with respect to x
double precision, intent(out),dimension(20)       :: dNidy !< Derivations of the shapefunctions with respect to y
double precision, intent(out),dimension(20)       :: dNidz !< Derivations of the shapefunctions with respect to z
double precision, intent(out)                     :: detjac !< Determinant of the jacobian matrix
double precision, intent(out)                     :: jac(3,3) !< Jacobi matrix
!
! Internal
!
double precision,dimension(3,3)               :: ijac
double precision,dimension(20)                :: dNidxi,dNideta,dNidzeta
integer                           :: j
!
integer                           :: err_code=0
logical                           :: ok_flag
!
! =================================================================================================
!
! Initialisation
!
  jac = 0.D0
!
! =================================================================================================
!
! Calculation
!
!  Compute derivations of the shape functions   
     call hexa20_shapefunctions_xi_eta_zeta(xi,eta,zeta,dNidxi=dNidxi,dNideta=dNideta,dNidzeta=dNidzeta)
!
!  Compute Jacobian
     do,j=1,20
       jac(1,1) = jac(1,1) + (dNidxi(j) * node_coordinates(j,1))
       jac(1,2) = jac(1,2) + (dNidxi(j) * node_coordinates(j,2))
       jac(1,3) = jac(1,3) + (dNidxi(j) * node_coordinates(j,3))
       jac(2,1) = jac(2,1) + (dNideta(j) * node_coordinates(j,1))
       jac(2,2) = jac(2,2) + (dNideta(j) * node_coordinates(j,2))
       jac(2,3) = jac(2,3) + (dNideta(j) * node_coordinates(j,3))
       jac(3,1) = jac(3,1) + (dNidzeta(j) * node_coordinates(j,1))
       jac(3,2) = jac(3,2) + (dNidzeta(j) * node_coordinates(j,2))
       jac(3,3) = jac(3,3) + (dNidzeta(j) * node_coordinates(j,3))
     end do
!       
!  Compute jacobian determinant (here only used for checking the element geometry)
!     detjac = det3x3(jac)  
!
!  Compute inverse jacobian (if jacobian determinant is positive)
!       if(detjac<1.d-12) then
! !   In case of wrong numbering or degenerated elements
!         write(*,*)
!         write(*,*) "Jacobi determinant is too small or negative! The geometry of one or more elements is invalid!"
!         err_code = 1
!         goto 9999
!       else
!        call inv3x3(jac,ijac)
        call M33INV (jac, ijac, detjac, OK_FLAG)
!       end if             
!
!  Compute derivations of the shapefunctions with respect to the global coordinates x,y,z
     do,j=1,20
       dNidx(j) = ijac(1,1) * dNidxi(j) + ijac(1,2) * dNideta(j) + ijac(1,3) * dNidzeta(j)
       dNidy(j) = ijac(2,1) * dNidxi(j) + ijac(2,2) * dNideta(j) + ijac(2,3) * dNidzeta(j)
       dNidz(j) = ijac(3,1) * dNidxi(j) + ijac(3,2) * dNideta(j) + ijac(3,3) * dNidzeta(j)
     end do
   
! =================================================================================================
!
! Error handling
!
9999 continue
!
if (err_code /= 0) then  
   write(*,*)                      'An error occured in subroutine'
   write(*,*)              'hexa20_shapefunctions_x_y_z'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
end if
return
end subroutine hexa20_shapefunctions_x_y_z
