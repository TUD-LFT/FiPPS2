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
subroutine quad8_geostiff_stress (disp_q8,Amat,Bmat,quad8_node_id_coords,quad8_Kg)
! =================================================================================================
!
!	Header:		control subroutine for computing element geometric stiffness matrix for flat
!			quadratic eight-node serendipity shell/plate element and assemble it to the
!			total geometric stiffness matrix
!
!	Content:	subroutine calculates the element geometric stiffness matrix for flat
!			quadratic eight-node quadrilateral shell/plate element and assembles it to
!			the total geometric stiffness matrix
!
!	Input:		
!
!	Output:		
!
!	Internal:	
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
!	Revision:	
!
! =================================================================================================
!
! Use
!
use netz_variablen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
!
! =================================================================================================
!
! Data types
!
! Input
!
double precision, dimension(3,3)		:: Amat,Bmat
double precision, dimension(8,4)		:: quad8_node_id_coords
!
! Output
!
double precision, dimension(48,48)		:: quad8_Kg
!
! Input + Output
!
double precision, dimension(48), intent(inout)	:: disp_q8
!
! inner
!
double precision, dimension(:,:)		:: TRMat(48,48),quad8_node_coords_lo(8,4)
double precision, dimension(:),allocatable	:: xi_vec, eta_vec
double precision, dimension(:),allocatable	:: w_xi,w_eta
double precision, dimension(8)			:: Ni
double precision, dimension(8)			:: dNidx, dNidy
double precision				:: detJac, facDP
double precision, dimension(:,:)		:: Bu(3,8),Bv(3,8),Bmem(3,16),Bbnd(3,16)
double precision, dimension(3)			:: nvec
double precision, dimension(8,8)		:: kg, Matd1d1, Matd2d2, Matd1d2!, Matd2d1
integer						:: inplaneintppd	! inplane integration point per dimension in natural coordinates
integer						:: numgausspoints

integer						:: ii, jj, kk, facI1, facI2
integer						:: err_code=0

!
! =================================================================================================
!
! Initialisation
!
! reduced integration ?
  inplaneintppd = 2
! initislise partial geometric stiffness matrix with zeros
  kg(1:8,1:8) = 0.D0
!
! =================================================================================================
!
! Calculation
!

! get transformation matrix from local to global coordinate system

  call quad8_rotation (quad8_node_id_coords,'lg',TRMat,quad8_node_coords_lo)
 
! transform displacement vector to local coordinates

  disp_q8=matmul(transpose(TRMat),disp_q8)

  numgausspoints = inplaneintppd*inplaneintppd                                                     ! total number of Gauss Points
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))
  
  call gauss_integration_values (inplaneintppd, xi_vec, eta_vec, w_xi, w_eta)

! loop over Gauss-points
   
  do ii=1,numgausspoints
   
    ! calculate shape-function values with respect to natural coordinates        ! == delta3

    call quad8_ansatzfunction(xi_vec(ii), eta_vec(ii), quad8_node_coords_lo, .false., Ni, dNidx, dNidy, detJac)

    do jj=1,8
       
       ! fill operator-matrices
       
       Bu(1,jj) = dNidx(jj)
       Bu(2,jj) = 0.D0
       Bu(3,jj) = dNidy(jj)
       
       Bv(1,jj) = 0.D0
       Bv(2,jj) = dNidy(jj)
       Bv(3,jj) = dNidx(jj)
    
    end do
    
    ! build matrices for u+v- and bex+bey-dofs
    
    Bmem(1:3,1:8)  = Bu
    Bmem(1:3,9:16) = Bv
    
    Bbnd(1:3,1:8)  = -Bv
    Bbnd(1:3,9:16) =  Bu
    
    ! multiply coefficient matrices with material membrane and membrane-bending-coupling stiffness amtrix
    
    Bmem = matmul(Amat,Bmem)
    Bbnd = matmul(Bmat,Bbnd)

    ! Calculate membrane forces
    
    do jj=1,3
       
       nvec(jj) = 0.D0
          
       do kk=1,size(Bmem,2)
          
          if (kk <= size(Bmem,2)/2) then
             
             facI1 = 0
             facI2 = 1
             
          else
             
             facI1 = size(Bmem,2)/2
             facI2 = 2
             
          end if
          
          facI1 = (kk-1-facI1)*ndof+facI2
          
          nvec(jj) = nvec(jj)+Bmem(jj,kk)*disp_q8(facI1)+Bbnd(jj,kk)*disp_q8(facI1+3)
          
       end do
       
    end do
    
    ! build matrices of dof derivative products
    
    do jj=1,8
       do kk=1,8
    
          Matd1d1(jj,kk) = dNidx(jj)*dNidx(kk)
          Matd2d2(jj,kk) = dNidy(jj)*dNidy(kk)
          Matd1d2(jj,kk) = dNidx(jj)*dNidy(kk)
    
       end do
    end do
    
    ! factor for numerical integration with Gauss-quadrature
    
    facDP = w_xi(ii)*w_eta(ii)*detJac
    
    ! Compute partial matrix kg with integration over element surface
    
  !  kg = kg+(nvec(1)*Matd1d1 + nvec(2)*Matd2d2 + 0.5D0*nvec(3)*(Matd1d2+transpose(Matd1d2)))*facDP
    kg = kg+(nvec(1)*Matd1d1 + nvec(2)*Matd2d2 + 1.0D0*nvec(3)*(Matd1d2+transpose(Matd1d2)))*facDP
   
  end do

  deallocate(xi_vec, eta_vec, w_xi, w_eta)

  ! Blow up kg to Kg and resort for nodal dofs on the fly

  do ii=1,8	  ! loop over rows
    do jj=1,8	  ! loop over columns
      do kk=1,3   ! loop over dofs u,v,w
      
	 quad8_Kg((ii-1)*ndof+kk,(jj-1)*ndof+kk) = quad8_Kg((ii-1)*ndof+kk,(jj-1)*ndof+kk)+kg(ii,jj)
	 
      end do
    end do
  end do

  ! Transform from local to global coordinates

  quad8_Kg = matmul(matmul(TRMat,quad8_Kg),transpose(TRMat))
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)  			   'quad8_geostiff_stress'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_geostiff_stress
