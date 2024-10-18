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
!> Computes the strain-displacement matrix for layered 20-node solid element
!> at a given point (xi,eta,zeta)
!
!> @details
!> Subroutine for computing the strain-displacement matrix for layered 20-node
!> solid element at a given point (xi,ety,zeta) with respect to global coordinate
!> system;
!> Ordering of the matrix corresponds to strain vector
!> e = [e1,e2,e3,gamma23,gamm13,gamma12]^T
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 29.06.2010
!
!> $Id: lsolid20_strain_displ_matrix.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_strain_displ_matrix(xi,eta,zeta,node_coords,Bmat,J,detjac)
! =================================================================================================
! use
!
! use konstanten
! use netz_variablen
! use globale_variablen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
! Input
!
double precision, intent(in)                   :: xi !< Natural coordinate xi at which strain-displacement matrix is calculated
double precision, intent(in)                   :: eta !< Natural coordinate eta at which strain-displacement matrix is calculated
double precision, intent(in)                   :: zeta !< Natural coordinate zeta at which strain-displacement matrix is calculated
double precision, intent(in)                   :: node_coords(20,3) !< Coordinates of element nodes in global coordinate system; dimension: (20,3)
!
! Output
double precision, intent(out), dimension(6,60) :: Bmat !< Strain-displacement matrix, as described above; dimension: (6,60)
double precision, intent(out), optional        :: detjac !< Determinant of Jacobi-Matrix (optional output)
double precision, intent(out), optional        :: J(3,3) !< Jacobi-Matrix (optional output); dimension: (3,3)
!
! Internal
!
double precision, dimension(20)                :: dNidx, dNidy, dNidz
double precision                               :: detJ, jac(3,3)
!
integer                                        :: err_code=0, ii
!
! =================================================================================================
!
! Calculation
!
 ! Get derivates of shape functions with respect to global coordinates
 call hexa20_shapefunctions_x_y_z(xi,eta,zeta,node_coords,dNidx,dNidy,dNidz,detJ,jac)

 ! Set up strain-displacement matrix
 Bmat = 0.d0
 
 do ii = 1,20
   Bmat(1,3*(ii-1)+1) = dNidx(ii)
   Bmat(2,3*(ii-1)+2) = dNidy(ii)
   Bmat(3,3*(ii-1)+3) = dNidz(ii)
   Bmat(4,3*(ii-1)+2) = dNidz(ii)
   Bmat(4,3*(ii-1)+3) = dNidy(ii)
   Bmat(5,3*(ii-1)+1) = dNidz(ii)
   Bmat(5,3*(ii-1)+3) = dNidx(ii)
   Bmat(6,3*(ii-1)+1) = dNidy(ii)
   Bmat(6,3*(ii-1)+2) = dNidx(ii)
 end do
 
 if (present(detjac)) then
   detjac = detJ
 end if
 
 if (present(J)) then
   J(:,:) = jac(:,:)
 end if

!
! ==========================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then

  write(*,*)                      'An error occured in subroutine'
  write(*,*)                      'lsolid20_strain_displ_matrix'
  write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
  write(*,*)                      'exit program '
  stop

end if

return

end subroutine lsolid20_strain_displ_matrix
