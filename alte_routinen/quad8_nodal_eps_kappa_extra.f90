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
subroutine quad8_nodal_eps_kappa_extra (disp_quad8, node_coords_lo, epsilon, kappa)

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
  integer                                         :: ii, ll
  integer                                         :: inplaneintppd, numgausspoints
  
  double precision				  :: detJac
  double precision, dimension(:), allocatable     :: xi_vec, eta_vec
  double precision, dimension(:),allocatable	  :: w_xi,w_eta
  double precision, dimension(8)                  :: Ni
  double precision, dimension(8)                  :: dNidx, dNidy
  double precision, dimension(3,8)                :: Bu, Bv
  double precision, dimension(:,:), allocatable   :: Transform
  double precision, dimension(:,:),allocatable	  :: preepsilon, prekappa

  
  integer                                         :: err_code = 0
  
  inplaneintppd = 2

  numgausspoints = inplaneintppd*inplaneintppd
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))

  allocate(Transform(8,numgausspoints))
  allocate(preepsilon(numgausspoints,3), prekappa(numgausspoints,3))
    
  call quad8_result_extrapolation_transform (inplaneintppd, Transform)
  
  call gauss_integration_values (inplaneintppd, xi_vec, eta_vec, w_xi, w_eta)
  
  epsilon = 0.d0
  kappa   = 0.d0
  preepsilon = 0.d0
  prekappa   = 0.d0

  do ii=1,numgausspoints ! loop over Gauss-points
  
    ! calculate shape-function values with respect to natural coordinates

    call quad8_ansatzfunction(xi_vec(ii), eta_vec(ii), node_coords_lo, .false., Ni, dNidx, dNidy, detJac)

    ! Calculate shape function derivatives with respect to local cartesian coordinate system
  
    ! fill operator-matrices
    
    Bu(1,1:8) = dNidx(1:8)
    Bu(2,1:8) = 0.d0
    Bu(3,1:8) = dNidy(1:8)
    
    Bv(1,1:8) = 0.D0
    Bv(2,1:8) = dNidy(1:8)
    Bv(3,1:8) = dNidx(1:8)
    
    do ll = 1,8
      preepsilon(ii,1:3) = preepsilon(ii,1:3) + ( Bu(:,ll) * disp_quad8((ll-1)*6 + 1) + Bv(:,ll) * disp_quad8((ll-1)*6 + 2) )
      prekappa(ii,1:3)   = prekappa(ii,1:3)   + ( Bu(:,ll) * disp_quad8((ll-1)*6 + 5) - Bv(:,ll) * disp_quad8((ll-1)*6 + 4) )
    end do
  
  end do ! End loop over Gauss-points

  epsilon(1:8,1:3) = matmul(Transform,preepsilon(1:numgausspoints,1:3))
  kappa(1:8,1:3)   = matmul(Transform,prekappa(1:numgausspoints,1:3))
  
  deallocate(xi_vec, eta_vec)
  deallocate(  w_xi,   w_eta)

  deallocate(Transform)
  deallocate(preepsilon, prekappa)
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_nodal_eps_kappa_extra'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return
  
end subroutine quad8_nodal_eps_kappa_extra
