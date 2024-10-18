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
subroutine quad8_iso_strains_gauss (elem, disp_quad8, pos, inplaneintppd, epsilon)
  
  use netz_variablen
  use globale_variablen, ONLY : mat1s

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
!
! Input
!
  integer, intent(in)                             :: pos
  double precision, dimension(:)                  :: disp_quad8(48)
  integer                                         :: inplaneintppd      ! Integrationspunkte in der Elementebene
!
! Output
! 
  double precision, dimension(:,:)                :: epsilon(inplaneintppd*inplaneintppd,6)
!
! Input + Output
!

!
! inner
!  
  integer                                         :: ii,jj,kk
  integer                                         :: elem
  double precision, dimension(:), allocatable     :: xi_vec, eta_vec, w_xi, w_eta
  double precision, dimension(:)                  :: reorderedDispl(48)
  double precision                                :: zeta
  double precision, dimension(8,3)                :: node_coords
  
  double precision, dimension(8)                  :: xi_n, eta_n
  
  integer                                         :: numgausspoints
  
  double precision, dimension(8)                  :: dNidxi,dNideta
  double precision                                :: dxdxi, dxdeta, dydxi, dydeta, dzdxi, dzdeta
  
  double precision, dimension(3,8)                :: normalVec              ! Normalenvektoren an den acht Punkten
  double precision, dimension(3)                  :: x_loc, y_loc, z_loc    ! lokale Koordinatenachsen
  double precision, dimension(8,3,3)              :: nodeTransMat           ! lokales Koordinatensystem an jedem Knoten
  double precision, dimension(3,3,8)              :: nodeTransMat2
  
  double precision, dimension(48,6)               :: bMat
  
  double precision                                :: detJac, dTemp
    
  integer                                         :: err_code = 0
  double precision, dimension(3,3)                :: transMat

!  double precision, dimension(3,3)                :: epsTensor
  double precision                                :: k
  
  xi_n(1) = -1.d0; eta_n(1) = -1.d0
  xi_n(2) =  1.d0; eta_n(2) = -1.d0
  xi_n(3) =  1.d0; eta_n(3) =  1.d0
  xi_n(4) = -1.d0; eta_n(4) =  1.d0
  xi_n(5) =  0.d0; eta_n(5) = -1.d0
  xi_n(6) =  1.d0; eta_n(6) =  0.d0
  xi_n(7) =  0.d0; eta_n(7) =  1.d0
  xi_n(8) = -1.d0; eta_n(8) =  0.d0
  
  numgausspoints = inplaneintppd*inplaneintppd
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))
  
  call gauss_integration_values (inplaneintppd, xi_vec, eta_vec, w_xi, w_eta)
  
  numgausspoints = inplaneintppd*inplaneintppd

  do kk = 1,8
    node_coords(kk,1:3) = nodes(quad8s(elem)%nids(kk))%coords(1:3)
  end do
  
  do jj = 1,6
    do ii = 1,8
      reorderedDispl((jj-1)*8+ii) = disp_quad8((ii-1)*6+jj)
    end do
  end do
  
  !*****************************************************
  ! Berechnung des Normalenvektors in den Knotenpunkten
  !*****************************************************
  do ii = 1,8
  
    call quad8_ansatzfunction_xieta(xi_n(ii), eta_n(ii), dNidxi=dNidxi, dNideta=dNideta)
    
    dxdxi  = 0.D0
    dxdeta = 0.D0
    dydxi  = 0.D0
    dydeta = 0.D0
    dzdxi  = 0.D0
    dzdeta = 0.D0
    
    do jj=1,8
    
      dxdxi  = dxdxi +dNidxi(jj) *node_coords(jj,1)
      dxdeta = dxdeta+dNideta(jj)*node_coords(jj,1)
    
      dydxi  = dydxi +dNidxi(jj) *node_coords(jj,2)
      dydeta = dydeta+dNideta(jj)*node_coords(jj,2)
    
      dzdxi  = dzdxi +dNidxi(jj) *node_coords(jj,3)
      dzdeta = dzdeta+dNideta(jj)*node_coords(jj,3)
    
    end do
    
    ! Berechnung der Normalen auf die Elementmittelsfläche
    normalVec(1,ii) = dydxi*dzdeta - dzdxi*dydeta
    normalVec(2,ii) = dzdxi*dxdeta - dxdxi*dzdeta
    normalVec(3,ii) = dxdxi*dydeta - dydxi*dxdeta
    
    normalVec(:,ii) = normalVec(:,ii) / SQRT(DOT_PRODUCT(normalVec(:,ii),normalVec(:,ii)))
        
    ! lokales Koordinatensystem berechnen
    
    x_loc(1) = dxdxi
    x_loc(2) = dydxi
    x_loc(3) = dzdxi
    
    x_loc = x_loc / SQRT(DOT_PRODUCT(x_loc,x_loc))
    
    ! Berechnung der Normalen auf die Elementmittelsfläche am Integrationspunkt
    z_loc(1:3) = normalVec(1:3,ii)
    
    ! cross_prod(z_vec, x_vec) (sollte normiert sein, da x und z normiert wurden)
    y_loc(1) = z_loc(2)*x_loc(3) - x_loc(2)*z_loc(3)
    y_loc(2) = z_loc(3)*x_loc(1) - x_loc(3)*z_loc(1)
    y_loc(3) = z_loc(1)*x_loc(2) - x_loc(1)*z_loc(2)

    nodeTransMat(ii,1:3,1) = -y_loc(1:3)
    nodeTransMat(ii,1:3,2) =  x_loc(1:3)
    nodeTransMat(ii,1:3,3) = 0.d0 
    
    nodeTransMat2(1:3,1,ii) = x_loc(1:3)
    nodeTransMat2(1:3,2,ii) = y_loc(1:3)
    nodeTransMat2(1:3,3,ii) = z_loc(1:3)
    
    nodeTransMat(ii,1:3,1:3) = matmul(nodeTransMat(ii,1:3,1:3),transpose(nodeTransMat2(1:3,1:3,ii)))
    
  end do

  zeta = 0.d0 ! Mitte
  if (pos .eq. 1) then !Oberseite
    zeta =  1.d0
  else if (pos .eq. 2) then !Unterseite
    zeta = -1.d0
  end if
  
  epsilon = 0.d0
  
  dTemp = mat1s(pshells(quad8s(elem)%int_pid)%intMat1ID)%nu
  dTemp = dTemp/(1.d0 - dTemp)
  
  k = MAX(1.2d0, 1.d0+0.2d0*quad8s(elem)%area/(25.d0*pshells(quad8s(elem)%int_pid)%mt**2))
  
  do ii=1,numgausspoints
    
    call quad8_iso_bmat(xi_vec(ii), eta_vec(ii), zeta, node_coords, pshells(quad8s(elem)%int_pid)%mt, normalVec, nodeTransMat, Bmat, detJac, transMat)

    epsilon(ii,1:6) = MatMul(transpose(Bmat),reorderedDispl)
    
    epsilon(ii,3)   = -dTemp * (epsilon(ii,1) + epsilon(ii,2))
    
    epsilon(ii,5)   = epsilon(ii,5)/k
    epsilon(ii,6)   = epsilon(ii,6)/k
    
!     epsTensor(1,1) =  epsilon(ii,1); epsTensor(1,2) = epsilon(ii,4)/2.d0; epsTensor(1,3) = epsilon(ii,6)/2.d0;
!     epsTensor(2,1) = epsTensor(1,2); epsTensor(2,2) = epsilon(ii,2)     ; epsTensor(2,3) = epsilon(ii,5)/2.d0;
!     epsTensor(3,1) = epsTensor(1,3); epsTensor(3,2) = epsTensor(2,3)    ; epsTensor(3,3) = epsilon(ii,3);
    
    !Transformation ins globale System
!    epsTensor = MATMUL(transMat,MATMUL(epsTensor,TRANSPOSE(transMat)))

    !WRITE(*,'(I8,X,6E12.5)') elem, 

  end do
  
  deallocate(xi_vec, eta_vec)
  deallocate(  w_xi,   w_eta)
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'quad8_iso_strains_gauss'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
    
  end if
  
  return
  
end subroutine quad8_iso_strains_gauss
