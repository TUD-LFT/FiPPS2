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
subroutine quad8_iso_results_node (elem, disp_quad8, pos, epsilon, sigma)
  
  use netz_variablen
  use globale_variablen, ONLY : mat1s

  implicit none
!
! =================================================================================================
!
! Input
!
  integer, intent(in)                             :: elem, pos
  double precision, dimension(:)		  :: disp_quad8(48)
!
! Output
! 
  double precision, dimension(:,:)                :: epsilon(8,6)
  double precision, dimension(:,:)                :: sigma(8,6)
!
! Input + Output
!

!
! inner
!  
  
  integer, parameter                              :: inplaneintppd  = 2		! Integrationspunkte in der Elementebene
  integer                                         :: numgausspoints
  
  double precision, dimension(:,:),allocatable	  :: preepsilon
  double precision, dimension(:,:),allocatable	  :: Transform
  
  integer                                         :: ii
  
  integer                                         :: err_code = 0
  
  double precision				  :: E, nue, thickness, fac
  double precision, dimension(6,6)                :: Dmat

  numgausspoints = inplaneintppd*inplaneintppd

  allocate(Transform(8,numgausspoints))
  allocate(preepsilon(numgausspoints,6))

  call quad8_iso_strains_gauss (elem, disp_quad8, pos, inplaneintppd, preepsilon)
  
  call quad8_result_extrapolation_transform (inplaneintppd, Transform)
  
  epsilon(1:8,1:6) = matmul(Transform,preepsilon(1:numgausspoints,1:6))
  
  ! Berechnen der D-Matrix
  
  E	    = mat1s(pshells(quad8s(elem)%int_pid)%intMat1ID)%ym
  nue	    = mat1s(pshells(quad8s(elem)%int_pid)%intMat1ID)%nu
  thickness = pshells(quad8s(elem)%int_pid)%mt
  
  fac = E / (1.d0 - nue**2)
  
  Dmat = 0.d0
  
  ! Hinweis: Der Schubkorrekturfaktor ist bereits auf den Dehnungen drauf.
  
  Dmat(1,1) =	  fac; Dmat(1,2) = nue*fac;
  Dmat(2,1) = nue*fac; Dmat(2,2) =     fac;
  Dmat(4,4) = fac * (1.d0 - nue)/2.d0
  Dmat(5,5) = fac * (1.d0 - nue)/2.d0
  Dmat(6,6) = fac * (1.d0 - nue)/2.d0
      
  DO ii = 1,8
  
    sigma(ii,1:6)    = MATMUL(Dmat,epsilon(ii,1:6))
        
  END DO
  
  deallocate(Transform)
  deallocate(preepsilon)
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_iso_strains_node'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return

end subroutine quad8_iso_results_node
