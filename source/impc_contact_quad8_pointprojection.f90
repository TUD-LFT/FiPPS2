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
!
!> @details
!
!> @author 
!
!> $Id: impc_contact_quad8_pointprojection.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine impc_contact_quad8_pointprojection(node_coords, xk, xi, eta)

  implicit none
!
! =================================================================================================
!
! Interface
!
  interface 
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
      double precision, intent(in)                              :: xi, eta
      
      double precision, dimension(8), optional, intent(out)     :: Ni
      double precision, dimension(8), optional, intent(out)     :: dNidxi, dNideta
      double precision, dimension(8), optional, intent(out)     :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
  end interface
  
  double precision, intent(out)    :: xi, eta
  
  double precision, dimension(8,3) :: node_coords
  double precision, dimension(2)   :: xieta
  double precision, dimension(3)   :: xk, xq, a1, a2
  double precision, dimension(3)   :: da1dxi, da1deta, da2dxi, da2deta
  double precision, dimension(8)   :: Ni
  double precision, dimension(8)   :: dNidxi, dNideta
  double precision, dimension(8)   :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
  double precision, dimension(2)   :: f
  double precision, dimension(2,2) :: Jac, invJac
  double precision                 :: detJac
  double precision, parameter      :: eps = 1.d-8
  integer, parameter               :: maxIter = 100
  
  integer                          :: ii, iter
  
  ! initial guess
  
  xieta = 0.0d0
  
  do iter = 1, maxIter
  
    call quad8_ansatzfunction_xieta(xieta(1), xieta(2), Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
    
    xq = 0.d0
    do ii = 1, 8
      xq(1:3) = xq(1:3) + Ni(ii)*node_coords(ii,1:3)
    end do
    
    a1 = 0.d0
    do ii = 1, 8
      a1(1:3) = a1(1:3) + dNidxi(ii)*node_coords(ii,1:3)
    end do
    
    a2 = 0.d0
    do ii = 1, 8
      a2(1:3) = a2(1:3) + dNideta(ii)*node_coords(ii,1:3)
    end do
    
    f(1) = SUM(a1*(xk-xq))
    f(2) = SUM(a2*(xk-xq))
    
    if (ABS(f(1)) .LT. eps .AND. ABS(f(2)) .LT. eps) exit
    
    da1dxi = 0.d0
    do ii = 1, 8
      da1dxi(1:3) = da1dxi(1:3) + d2Nidxi2(ii)*node_coords(ii,1:3)
    end do
    
    da1deta = 0.d0
    do ii = 1, 8
      da1deta(1:3) = da1deta(1:3) + d2Nidxideta(ii)*node_coords(ii,1:3)
    end do
    
    da2dxi = 0.d0
    do ii = 1, 8
      da2dxi(1:3) = da2dxi(1:3) + d2Nidetadxi(ii)*node_coords(ii,1:3)
    end do
    
    da2deta = 0.d0
    do ii = 1, 8
      da2deta(1:3) = da2deta(1:3) + d2Nideta2(ii)*node_coords(ii,1:3)
    end do
    
    Jac(1,1) = SUM(xk*da1dxi)  - SUM(xq*da1dxi)  - SUM(a1**2) ! dF1dxi
    Jac(1,2) = SUM(xk*da1deta) - SUM(xq*da1deta) - SUM(a1*a2) ! dF1deta
    Jac(2,1) = SUM(xk*da2dxi)  - SUM(xq*da2dxi)  - SUM(a1*a2) ! dF2dxi
    Jac(2,2) = SUM(xk*da2deta) - SUM(xq*da2deta) - SUM(a2**2) ! dF2deta
    
    detJac = 1.d0/(Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1))
    
    invJac(1,1) =  detJac*Jac(2,2)
    invJac(1,2) = -detJac*Jac(1,2)
    invJac(2,1) = -detJac*Jac(2,1)
    invJac(2,2) =  detJac*Jac(1,1)
    
    xieta = xieta - MatMul(invJac,f)
    
  end do
  
  if (iter .eq. maxIter) then
    WRITE(*,*) 'Kein Projektionspunkt innerhalb der max. Iterationsanzahl gefunden!'
  else if (ABS(xieta(1)) .GT. 1.0000001D0 .OR. ABS(xieta(2)) .GT. 1.0000001d0) then
    WRITE(*,*) 'Der gefundene Projektionspunkt liegt ausserhalb des Elements!', xieta(1), xieta(2)
  end if
  
  xi  = MAX(MIN(xieta(1),1.d0), -1.d0)
  eta = MAX(MIN(xieta(2),1.d0), -1.d0)
  
end subroutine impc_contact_quad8_pointprojection
