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
!  Get natural coordinates of lsolid20-elements for a specific point
!
!> @details
!  The natural coordinates of the layered solid 20-node element defined by the nodes given
!  in node_coords are calculated for the point xk given in global coordinates. Therefore,
!  the least squares problem minimizing the distance between xk and the point described by
!  the resulting natural coordinates xi, eta, zeta is solved numerically using Newtons method.
!
!> @author
!> Florian Dexl, TU Dresden, wiss. Mitarbeiter
!
!> @date
!  21.04.2022
!
!> $Id: impc_contact_quad8_pointprojection.f90 413 2019-10-17 08:20:58Z DOM+ahauffe $
!> $Author: DOM+ahauffe $
!> $Revision: 413 $
!> $Date: 2019-10-17 10:20:58 +0200 (Do, 17. Okt 2019) $
!
! =================================================================================================
subroutine impc_contact_lsolid20_naturalcoordinates(node_coords, xk, xi, eta, zeta)
! use
!
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
  
  double precision, intent(out)     :: xi, eta, zeta
  
  double precision, dimension(20,3) :: node_coords
  double precision, dimension(3)    :: xietazeta
  double precision, dimension(3)    :: xk, xq
  double precision, dimension(20)   :: Ni
  double precision, dimension(20)   :: dNidxi, dNideta, dNidzeta
  double precision, dimension(3)    :: f
  double precision, dimension(3,3)  :: Jac, invJac
  double precision                  :: detJac
  double precision, parameter       :: eps = 1.d-8
  integer, parameter                :: maxIter = 100
  
  integer                           :: ii, iter

  integer                           :: err_code=0
  logical                           :: ok_flag
  
  ! initial guess
  xietazeta(:) = 0.0d0
  
  do iter = 1, maxIter
  
    ! Get shape functions and first derivatives
    call hexa20_shapefunctions_xi_eta_zeta(xietazeta(1), xietazeta(2), xietazeta(3), Ni, dNidxi, dNideta, dNidzeta)

    ! Get global coordinates for current iteration
    xq(:) = 0.d0
    do ii = 1, 20
      xq(:) = xq(:) + Ni(ii)*node_coords(ii,1:3)
    end do

    ! Objective is to solve the least squares problem minimizing the
    ! difference between xk(:) and xq(:) --> f is the residual function
    ! being minimized
    f(:) = xk(:) - xq(:)

    ! Check for convergence
    if (sqrt(sum(f(:)**2)) .lt. eps) exit

    ! Build Jacobi-matrix
    Jac(:,:) = 0.d0
    do ii = 1, 20
        ! First column: derivative of residual function f with respect to natural coordinate xi
        Jac(:,1) = Jac(:,1) - dNidxi(ii)*node_coords(ii,1:3)
        ! Second column: derivative of residual function f with respect to natural coordinate eta
        Jac(:,2) = Jac(:,2) - dNideta(ii)*node_coords(ii,1:3)
        ! Third column: derivative of residual function f with respect to natural coordinate zeta
        Jac(:,3) = Jac(:,3) - dNidzeta(ii)*node_coords(ii,1:3)
    end do

    ! Get inverse of Jacobi-matrix
    call M33INV (Jac, invJac, detJac, ok_flag)

    ! Check for error during matrix inversion
    if (ok_flag .eq. .false.) goto 9999

    ! Calulate new natural coordinates
    xietazeta(:) = xietazeta(:) - matmul(invJac,f)
    
  end do
  
  if (iter .eq. maxIter) then
    WRITE(*,*) 'Keine natuerlichen Koordinaten innerhalb der max. Iterationsanzahl gefunden!'
  else if (ABS(xietazeta(1)) .GT. 1.0000001D0 .OR. ABS(xietazeta(2)) .GT. 1.0000001d0 .OR. ABS(xietazeta(3)) .GT. 1.0000001d0) then
    WRITE(*,*) 'Der zu kontaktierende Knoten liegt ausserhalb des Lsolid20-Elements!', xietazeta(1), xietazeta(2), xietazeta(3)
  end if
  
  xi   = xietazeta(1)
  eta  = xietazeta(2)
  zeta = xietazeta(3)
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'impc_contact_lsolid20_naturalcoordinates'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return
  
end subroutine impc_contact_lsolid20_naturalcoordinates
