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
!> Subroutine for computing the shapefunctions and their first order
!> derivations with respect to the natural coordinates xi, eta, zeta.
!
!> @details
!> Computes the shapefunctions and the derivations with respect to 
!> the natural element coordinate system. Node numbering is equal to ANSYS
!> SOLID95.
!
!> @author Thomas Dylla, TU Dresden, Diplomarbeit, 13.12.2013
!
!> $Id: hexa20_shapefunctions_xi_eta_zeta.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine hexa20_shapefunctions_xi_eta_zeta(xi,eta,zeta,Ni,dNidxi,dNideta,dNidzeta)
! =================================================================================================
! use
!
use globale_variablen
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
double precision, intent(in)                  :: xi !< Natural element coordinate xi
double precision, intent(in)                  :: eta !< Natural element coordinate eta
double precision, intent(in)                  :: zeta !< Natural element coordinate zeta
!
! Output
!
double precision, optional ,intent(out),dimension(20) :: Ni !< Shapefunctions of the element at xi, eta, zeta
double precision, optional ,intent(out),dimension(20) :: dNidxi !< First derivative of the shapefunctions with respect to xi
double precision, optional ,intent(out),dimension(20) :: dNideta !< First derivative of the shapefunctions with respect to eta
double precision, optional ,intent(out),dimension(20) :: dNidzeta !< First derivative of the shapefunctions with respect to zeta
!
! Internal
!
integer                           :: err_code=0
!
! =================================================================================================
!
! Calculation
!
!    Compute shapefunctions (optional)
      if(present(Ni) .eq. .true.) then
    Ni(1) = .125D0 * (1.D0-xi) * (1.D0-eta) * (1.D0-zeta) * (-xi-eta-zeta-2.D0)
    Ni(2) = .125D0 * (1.D0+xi) * (1.D0-eta) * (1.D0-zeta) * ( xi-eta-zeta-2.D0)
    Ni(3) = .125D0 * (1.D0+xi) * (1.D0+eta) * (1.D0-zeta) * ( xi+eta-zeta-2.D0)
    Ni(4) = .125D0 * (1.D0-xi) * (1.D0+eta) * (1.D0-zeta) * (-xi+eta-zeta-2.D0)
    Ni(5) = .125D0 * (1.D0-xi) * (1.D0-eta) * (1.D0+zeta) * (-xi-eta+zeta-2.D0)
    Ni(6) = .125D0 * (1.D0+xi) * (1.D0-eta) * (1.D0+zeta) * ( xi-eta+zeta-2.D0)
    Ni(7) = .125D0 * (1.D0+xi) * (1.D0+eta) * (1.D0+zeta) * ( xi+eta+zeta-2.D0)
    Ni(8) = .125D0 * (1.D0-xi) * (1.D0+eta) * (1.D0+zeta) * (-xi+eta+zeta-2.D0)
!   
    Ni(9)  = .25D0 * (1.D0-xi*xi)*(1.D0-eta)*(1.D0-zeta)
    Ni(11) = .25D0 * (1.D0-xi*xi)*(1.D0+eta)*(1.D0-zeta)
    Ni(13) = .25D0 * (1.D0-xi*xi)*(1.D0-eta)*(1.D0+zeta)
    Ni(15) = .25D0 * (1.D0-xi*xi)*(1.D0+eta)*(1.D0+zeta)
!   
    Ni(10) = .25D0 * (1.D0-eta*eta)*(1.D0+xi)*(1.D0-zeta)
    Ni(12) = .25D0 * (1.D0-eta*eta)*(1.D0-xi)*(1.D0-zeta)
    Ni(14) = .25D0 * (1.D0-eta*eta)*(1.D0+xi)*(1.D0+zeta)
    Ni(16) = .25D0 * (1.D0-eta*eta)*(1.D0-xi)*(1.D0+zeta)
!   
    Ni(17) = .25D0 * (1.D0-zeta*zeta)*(1.D0-eta)*(1.D0-xi)
    Ni(18) = .25D0 * (1.D0-zeta*zeta)*(1.D0-eta)*(1.D0+xi)
    Ni(19) = .25D0 * (1.D0-zeta*zeta)*(1.D0+eta)*(1.D0+xi)
    Ni(20) = .25D0 * (1.D0-zeta*zeta)*(1.D0+eta)*(1.D0-xi)
      end if
!
!    Compute derivations of the shapefunctions with respect to xi (optional)
      if(present(dNidxi) .eq. .true.) then
         dNidxi(1) = .125D0 * (1.D0-eta)*(1.D0-zeta)*(-(-2.D0*xi-eta-zeta-1.D0))
         dNidxi(2) = .125D0 * (1.D0-eta)*(1.D0-zeta)*( ( 2.D0*xi-eta-zeta-1.D0))
         dNidxi(3) = .125D0 * (1.D0+eta)*(1.D0-zeta)*( ( 2.D0*xi+eta-zeta-1.D0))
         dNidxi(4) = .125D0 * (1.D0+eta)*(1.D0-zeta)*(-(-2.D0*xi+eta-zeta-1.D0))
      
     dNidxi(5) = .125D0 * (1.D0-eta)*(1.D0+zeta)*(-(-2.D0*xi-eta+zeta-1.D0))
         dNidxi(6) = .125D0 * (1.D0-eta)*(1.D0+zeta)*( ( 2.D0*xi-eta+zeta-1.D0))
         dNidxi(7) = .125D0 * (1.D0+eta)*(1.D0+zeta)*( ( 2.D0*xi+eta+zeta-1.D0))
         dNidxi(8) = .125D0 * (1.D0+eta)*(1.D0+zeta)*(-(-2.D0*xi+eta+zeta-1.D0))
      
         dNidxi(9 ) = -.5D0*xi*(1.D0-eta)*(1.D0-zeta)
         dNidxi(11) = -.5D0*xi*(1.D0+eta)*(1.D0-zeta)
         dNidxi(13) = -.5D0*xi*(1.D0-eta)*(1.D0+zeta)
         dNidxi(15) = -.5D0*xi*(1.D0+eta)*(1.D0+zeta)
         
         dNidxi(10) =  .25D0*(1.D0-eta*eta)*(1.D0-zeta)
         dNidxi(12) = -.25D0*(1.D0-eta*eta)*(1.D0-zeta)
         dNidxi(14) =  .25D0*(1.D0-eta*eta)*(1.D0+zeta)
         dNidxi(16) = -.25D0*(1.D0-eta*eta)*(1.D0+zeta)
      
     dNidxi(17) = -.25D0*(1.D0-zeta*zeta)*(1.D0-eta)
     dNidxi(18) =  .25D0*(1.D0-zeta*zeta)*(1.D0-eta)
     dNidxi(19) =  .25D0*(1.D0-zeta*zeta)*(1.D0+eta)
     dNidxi(20) = -.25D0*(1.D0-zeta*zeta)*(1.D0+eta) 
      end if
!   
!    Compute derivations of the shapefunctions with respect to eta (optional)
      if(present(dNideta) .eq. .true.) then
     dNideta(1) = .125D0 * (1.D0-xi)*(1.D0-zeta)*(-(-2.D0*eta-xi-zeta-1.D0))
         dNideta(2) = .125D0 * (1.D0+xi)*(1.D0-zeta)*(-(-2.D0*eta+xi-zeta-1.D0))
         dNideta(3) = .125D0 * (1.D0+xi)*(1.D0-zeta)*( ( 2.D0*eta+xi-zeta-1.D0))
         dNideta(4) = .125D0 * (1.D0-xi)*(1.D0-zeta)*( ( 2.D0*eta-xi-zeta-1.D0))
      
     dNideta(5) = .125D0 * (1.D0-xi)*(1.D0+zeta)*(-(-2.D0*eta-xi+zeta-1.D0))
         dNideta(6) = .125D0 * (1.D0+xi)*(1.D0+zeta)*(-(-2.D0*eta+xi+zeta-1.D0))
         dNideta(7) = .125D0 * (1.D0+xi)*(1.D0+zeta)*( ( 2.D0*eta+xi+zeta-1.D0))
         dNideta(8) = .125D0 * (1.D0-xi)*(1.D0+zeta)*( ( 2.D0*eta-xi+zeta-1.D0))
               
         dNideta(9 ) = -.25D0*(1.D0-xi*xi)*(1.D0-zeta)
         dNideta(11) =  .25D0*(1.D0-xi*xi)*(1.D0-zeta)
         dNideta(13) = -.25D0*(1.D0-xi*xi)*(1.D0+zeta)
         dNideta(15) =  .25D0*(1.D0-xi*xi)*(1.D0+zeta)
      
         dNideta(10) = -.5D0*eta*(1.D0+xi)*(1.D0-zeta)
         dNideta(12) = -.5D0*eta*(1.D0-xi)*(1.D0-zeta)
         dNideta(14) = -.5D0*eta*(1.D0+xi)*(1.D0+zeta)
         dNideta(16) = -.5D0*eta*(1.D0-xi)*(1.D0+zeta)

     dNideta(17) = -.25D0*(1.D0-zeta*zeta)*(1.D0-xi)
     dNideta(18) = -.25D0*(1.D0-zeta*zeta)*(1.D0+xi)
     dNideta(19) =  .25D0*(1.D0-zeta*zeta)*(1.D0+xi)
     dNideta(20) =  .25D0*(1.D0-zeta*zeta)*(1.D0-xi)    
      end if
!
!    Compute derivations of the shapefunctions with respect to zeta (optional)
      if(present(dNidzeta) .eq. .true.) then
         dNidzeta(1) = .125D0 * (1.D0-eta)*(1.D0-xi)*(-(-2.D0*zeta-eta-xi-1.D0))
         dNidzeta(2) = .125D0 * (1.D0-eta)*(1.D0+xi)*(-(-2.D0*zeta-eta+xi-1.D0))
         dNidzeta(3) = .125D0 * (1.D0+eta)*(1.D0+xi)*(-(-2.D0*zeta+eta+xi-1.D0))
         dNidzeta(4) = .125D0 * (1.D0+eta)*(1.D0-xi)*(-(-2.D0*zeta+eta-xi-1.D0))
      
     dNidzeta(5) = .125D0 * (1.D0-eta)*(1.D0-xi)*( ( 2.D0*zeta-eta-xi-1.D0))
         dNidzeta(6) = .125D0 * (1.D0-eta)*(1.D0+xi)*( ( 2.D0*zeta-eta+xi-1.D0))
         dNidzeta(7) = .125D0 * (1.D0+eta)*(1.D0+xi)*( ( 2.D0*zeta+eta+xi-1.D0))
         dNidzeta(8) = .125D0 * (1.D0+eta)*(1.D0-xi)*( ( 2.D0*zeta+eta-xi-1.D0))
       
     dNidzeta(9 ) = -.25D0*(1.D0-xi*xi)*(1.D0-eta)
     dNidzeta(11) = -.25D0*(1.D0-xi*xi)*(1.D0+eta)
     dNidzeta(13) =  .25D0*(1.D0-xi*xi)*(1.D0-eta)
     dNidzeta(15) =  .25D0*(1.D0-xi*xi)*(1.D0+eta)
       
         dNidzeta(10) = -.25D0*(1.D0-eta*eta)*(1.D0+xi)
         dNidzeta(12) = -.25D0*(1.D0-eta*eta)*(1.D0-xi)
         dNidzeta(14) =  .25D0*(1.D0-eta*eta)*(1.D0+xi)
         dNidzeta(16) =  .25D0*(1.D0-eta*eta)*(1.D0-xi)
      
         dNidzeta(17) = -.5D0*zeta*(1.D0-eta)*(1.D0-xi)
         dNidzeta(18) = -.5D0*zeta*(1.D0-eta)*(1.D0+xi)
         dNidzeta(19) = -.5D0*zeta*(1.D0+eta)*(1.D0+xi)
         dNidzeta(20) = -.5D0*zeta*(1.D0+eta)*(1.D0-xi)
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
   write(*,*)              'hexa20_shapefunctions_xi_eta_zeta'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
!   
end if
!
return
!
end subroutine hexa20_shapefunctions_xi_eta_zeta
