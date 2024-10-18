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
subroutine lsolid20_integration_values(nip, nop, xi_vec, eta_vec, zeta_vec, w_xi, w_eta, w_zeta)
! =================================================================================================
!
!   Header:     Calculate integration points and weights in 3 dimensions for layered 20-node
!               solid element
!
!   Content:    Returns the coordinates of the integration points and the corredponding weights
!               for the layered 20-node solid element; integration points in plane are given
!               by Gauss integration rule and 3, 2 or 1 point in each of the tewo directions
!               can be chosen; integration points through the thickness are given by Simpson's
!               integration rule and 9, 7, 5, 3 or 1 point in thickness direction can be chosen
!
!   Input:      nip:      number of integration points in plane (per direction) (3, 2 or 1)
!               nop:      number of integration points over thickness (9, 7, 5, 3 or 1)
!
!   Output:     xi_vec:   array containing xi-coordinates of integration points
!               eta_vec:  array containing eta-coordinates of integration points
!               zeta_vec: array containing zeta-coordinates of integration points
!               w_xi:     weights needed for integration over xi
!               w_eta:    weights needed for integration over eta
!               w_zeta:   weights needed for integration over zeta
!
!   Calls:      gauss_integration_values
!
!   Called by:  lsolid20_initstressvec
!               lsolid20_stiff
!               lsolid20_temploadvec
!               lsolid20_strains_stresses
!
!   Author:     Florian Dexl
!               TU Dresden, Diplomarbeit 2015
!
! =================================================================================================
!
! use
!
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
integer, intent(in)                                    :: nip, nop
!
! Output
!
double precision, dimension(nip*nip*nop), intent(out)  :: xi_vec, eta_vec, zeta_vec
double precision, dimension(nip*nip*nop), intent(out)  :: w_xi, w_eta, w_zeta
!
! Internal
!
integer                                                :: err_code=0, nip2, ii
!
! =================================================================================================
!
! Calculation
!
 ! Get in-plane integration points by Gauss rule
 nip2 = nip*nip
 call gauss_integration_values(nip, xi_vec(1:nip2), eta_vec(1:nip2), w_xi(1:nip2), w_eta(1:nip2))

 ! Get through thickness integration points by Simpson's rule
 if (nop == 1) then
   zeta_vec = 0.d0
   
   w_zeta   = 2.d0
 
 else if (nop == 3) then
 
   do ii = 1,(nop-1)
     xi_vec( ii*nip2+1:(ii+1)*nip2) = xi_vec(1:nip2)
     eta_vec(ii*nip2+1:(ii+1)*nip2) = eta_vec(1:nip2)
     
     w_xi(   ii*nip2+1:(ii+1)*nip2) = w_xi(1:nip2)
     w_eta(  ii*nip2+1:(ii+1)*nip2) = w_eta(1:nip2)
   end do
     
   zeta_vec(       1:       nip2) = -1.d0
   zeta_vec(    nip2+1:     2*nip2) =  0.d0
   zeta_vec(  2*nip2+1:     3*nip2) =  1.d0
   
   w_zeta(         1:       nip2) = 1.d0/3.d0
   w_zeta(      nip2+1:     2*nip2) = 4.d0/3.d0
   w_zeta(    2*nip2+1:     3*nip2) = 1.d0/3.d0 

 else if (nop == 5) then
 
   do ii = 1,(nop-1)
     xi_vec( ii*nip2+1:(ii+1)*nip2) = xi_vec(1:nip2)
     eta_vec(ii*nip2+1:(ii+1)*nip2) = eta_vec(1:nip2)
     
     w_xi(   ii*nip2+1:(ii+1)*nip2) = w_xi(1:nip2)
     w_eta(  ii*nip2+1:(ii+1)*nip2) = w_eta(1:nip2)
   end do

   zeta_vec(       1:       nip2) = -1.d0
   zeta_vec(    nip2+1:     2*nip2) = -0.5d0
   zeta_vec(  2*nip2+1:     3*nip2) =  0.d0
   zeta_vec(  3*nip2+1:     4*nip2) =  0.5d0
   zeta_vec(  4*nip2+1:     5*nip2) =  1.d0
   
   w_zeta(         1:       nip2) = 1.d0/6.d0
   w_zeta(      nip2+1:     2*nip2) = 4.d0/6.d0
   w_zeta(    2*nip2+1:     3*nip2) = 2.d0/6.d0
   w_zeta(    3*nip2+1:     4*nip2) = 4.d0/6.d0
   w_zeta(    4*nip2+1:     5*nip2) = 1.d0/6.d0

 else if (nop == 7) then

   do  ii = 1,(nop-1)
     xi_vec( ii*nip2+1:(ii+1)*nip2) = xi_vec(1:nip2)
     eta_vec(ii*nip2+1:(ii+1)*nip2) = eta_vec(1:nip2)
     
     w_xi(   ii*nip2+1:(ii+1)*nip2) = w_xi(1:nip2)
     w_eta(  ii*nip2+1:(ii+1)*nip2) = w_eta(1:nip2)
   end do

   zeta_vec(       1:       nip2) = -1.d0
   zeta_vec(    nip2+1:     2*nip2) = -2.d0/3.d0
   zeta_vec(  2*nip2+1:     3*nip2) = -1.d0/3.d0
   zeta_vec(  3*nip2+1:     4*nip2) =  0.d0
   zeta_vec(  4*nip2+1:     5*nip2) =  1.d0/3.d0
   zeta_vec(  5*nip2+1:     6*nip2) =  2.d0/3.d0
   zeta_vec(  6*nip2+1:     7*nip2) =  1.d0
   
   w_zeta(         1:       nip2) = 1.d0/9.d0
   w_zeta(      nip2+1:     2*nip2) = 4.d0/9.d0
   w_zeta(    2*nip2+1:     3*nip2) = 2.d0/9.d0
   w_zeta(    3*nip2+1:     4*nip2) = 4.d0/9.d0
   w_zeta(    4*nip2+1:     5*nip2) = 2.d0/9.d0
   w_zeta(    5*nip2+1:     6*nip2) = 4.d0/9.d0
   w_zeta(    6*nip2+1:     7*nip2) = 1.d0/9.d0

 else if (nop == 9) then

   do  ii = 1,(nop-1)
     xi_vec( ii*nip2+1:(ii+1)*nip2) = xi_vec(1:nip2)
     eta_vec(ii*nip2+1:(ii+1)*nip2) = eta_vec(1:nip2)
     
     w_xi(   ii*nip2+1:(ii+1)*nip2) = w_xi(1:nip2)
     w_eta(  ii*nip2+1:(ii+1)*nip2) = w_eta(1:nip2)
   end do

   zeta_vec(       1:       nip2) = -1.d0
   zeta_vec(    nip2+1:     2*nip2) = -3.d0/4.d0
   zeta_vec(  2*nip2+1:     3*nip2) = -2.d0/4.d0
   zeta_vec(  3*nip2+1:     4*nip2) = -1.d0/4.d0
   zeta_vec(  4*nip2+1:     5*nip2) =  0.d0
   zeta_vec(  5*nip2+1:     6*nip2) =  1.d0/4.d0
   zeta_vec(  6*nip2+1:     7*nip2) =  2.d0/4.d0
   zeta_vec(  7*nip2+1:     8*nip2) =  3.d0/4.d0
   zeta_vec(  8*nip2+1:     9*nip2) =  1.d0
   
   w_zeta(         1:       nip2) = 1.d0/12.d0
   w_zeta(      nip2+1:     2*nip2) = 4.d0/12.d0
   w_zeta(    2*nip2+1:     3*nip2) = 2.d0/12.d0
   w_zeta(    3*nip2+1:     4*nip2) = 4.d0/12.d0
   w_zeta(    4*nip2+1:     5*nip2) = 2.d0/12.d0
   w_zeta(    5*nip2+1:     6*nip2) = 4.d0/12.d0
   w_zeta(    6*nip2+1:     7*nip2) = 2.d0/12.d0
   w_zeta(    7*nip2+1:     8*nip2) = 4.d0/12.d0
   w_zeta(    8*nip2+1:     9*nip2) = 1.d0/12.d0
 
 else
   
   write(*,*) 'wrong input on parameter nop (number of'
   write(*,*) 'integration points over layer thickness)'
   write(*,*) 'values 1, 3, 5, 7 or 9 are allowed to be chosen'
   err_code = 2
   goto 9999

 end if

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
   write(*,*)                      'lsolid20_integration_values'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!   
end if
!
return
!
end subroutine lsolid20_integration_values