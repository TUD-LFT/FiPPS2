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
subroutine quad8_nodal_eps_kappa (disp_quad8, node_coords_lo, epsilon, kappa)

! =================================================================================================
!
!	Header:		
!
!	Content:	
!
!	Input:		
!
!	Output:		
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Andreas Hauffe			04.10.2010
! 			TU Dresden, WiMi
!
! =================================================================================================
!
use netz_variablen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Input
!
  double precision, dimension(8*ndof), intent(in) :: disp_quad8
  double precision, dimension(8,4), intent(in)    :: node_coords_lo
! Output
! 
  double precision, dimension(8,3), intent(out)   :: epsilon, kappa
!
! Input + Output
!

!
! inner
!
  integer                                         :: node
  integer                                         :: ll
  
  double precision				  :: detJac
  double precision, dimension(8)                  :: xi_vec, eta_vec
  double precision, dimension(8)		  :: Ni
  double precision, dimension(8)		  :: dNidx, dNidy
  double precision, dimension(3,8)                :: Bu, Bv
  
  integer                                         :: err_code = 0

    
  xi_vec(1) =  -1.d0; eta_vec(1) =  -1.d0
  xi_vec(2) =	1.d0; eta_vec(2) =  -1.d0
  xi_vec(3) =	1.d0; eta_vec(3) =   1.d0
  xi_vec(4) =  -1.d0; eta_vec(4) =   1.d0
  xi_vec(5) =  0.0d0; eta_vec(5) = -1.0d0
  xi_vec(6) =  1.0d0; eta_vec(6) =  0.0d0
  xi_vec(7) =  0.0d0; eta_vec(7) =  1.0d0
  xi_vec(8) = -1.0d0; eta_vec(8) =  0.0d0
  
  epsilon = 0.d0
  kappa   = 0.d0
  
  do node = 1,8
  
    call quad8_ansatzfunction(xi_vec(node), eta_vec(node), node_coords_lo, .false., Ni, dNidx, dNidy, detJac)
    
    Bu(1,1:8) = dNidx(1:8)
    Bu(2,1:8) = 0.d0
    Bu(3,1:8) = dNidy(1:8)
    
    Bv(1,1:8) = 0.D0
    Bv(2,1:8) = dNidy(1:8)
    Bv(3,1:8) = dNidx(1:8)
    
    do ll = 1,8
    
      epsilon(node,1:3) = epsilon(node,1:3) + Bu(:,ll) * disp_quad8((ll-1)*6 + 1) + Bv(:,ll) * disp_quad8((ll-1)*6 + 2)
      kappa(node,1:3)	= kappa(node,1:3)   + Bu(:,ll) * disp_quad8((ll-1)*6 + 5) - Bv(:,ll) * disp_quad8((ll-1)*6 + 4)
    
    end do
    
  end do
    
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_nodal_eps_kappa'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return
  
end subroutine quad8_nodal_eps_kappa
