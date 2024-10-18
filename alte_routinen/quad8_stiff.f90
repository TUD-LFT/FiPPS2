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
subroutine quad8_stiff (Amatel,Bmatel,Dmatel,Hmatel,SOmega,quad8_node_id_coords,quad8_K, area, thick)

! =================================================================================================
!
!	Header:		control subroutine for computing element stiffness matrix for linear 
! 			four-node quadrilateral shell/plate element and assemble it to the total
! 			stiffness matrix
!
!	Content:	subroutine calculates the element stiffness matrix for linear four-node
! 			quadrilateral shell/plate element and assembles it to the total stiffness matrix
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
use globale_variablen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Data types
!
! Input
!
double precision, dimension(3,3), intent(in)	:: Amatel, Bmatel, Dmatel
double precision, dimension(:,:), intent(in)	:: Hmatel(2,2), quad8_node_id_coords(8,4)
!
! Output
!
double precision, dimension(:,:), intent(out)	:: quad8_K(48,48)
double precision, intent(out)                   :: area
!
! Input + Output
!
double precision, intent(in)			:: SOmega
!
! inner
!
double precision, dimension(:,:)		:: TRMat(48,48), quad8_node_coords_lo(8,4)
double precision, dimension(:),allocatable	:: xi_vec, eta_vec
double precision, dimension(:),allocatable	:: w_xi,w_eta
double precision, dimension(8)			:: Ni
double precision, dimension(8)			:: dNidx, dNidy
double precision				:: detJac
double precision, dimension(:,:)		:: Bu(3,8),Bv(3,8),Bw(2,8),Bbex(2,8),Bbey(2,8),fac8t3(8,3),fac8t2(8,2)
double precision				:: SOmegaT
double precision                                :: thick

double precision, dimension(:,:)         	:: Hmatelt(2,2)

double precision, dimension(8,8)		:: kuu,kuv,kubex,kubey,kubez,kvv,kvbex,kvbey,kvbez
double precision, dimension(8,8)		:: kww,kwbex,kwbey,kbexbex,kbexbey,kbeybey,kbezbez

integer						:: ii,jj,kk
double precision				:: facDP
logical						:: int_selective
integer						:: inplaneintppd	! inplane integration point per dimension in natural coordinates
integer						:: inplaneintppdts
integer						:: numgausspoints

integer						:: err_code=0
!
! =================================================================================================
!
! Initialisation
!
  kuu=0.D0; kuv=0.D0; kubex=0.D0; kubey=0.D0; kubez=0.D0
  	    kvv=0.D0; kvbex=0.D0; kvbey=0.D0; kvbez=0.D0
  kww=0.D0; kwbex=0.D0; kwbey=0.D0
  kbexbex=0.D0; kbexbey=0.D0; kbeybey=0.D0; kbezbez=0.D0

  quad8_K = 0.D0

! selective reduced integration ?

  int_selective   = .false.
  inplaneintppd   = 2
  inplaneintppdts = 2

!
! =================================================================================================
!
! Calculation
!
! get transformation matrix from local to global coordinate system

  call quad8_rotation (quad8_node_id_coords,'lg',TRMat,quad8_node_coords_lo)
  
  numgausspoints = inplaneintppd*inplaneintppd                                                     ! total number of Gauss Points
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))
  
  call gauss_integration_values (inplaneintppd, xi_vec, eta_vec, w_xi, w_eta)
   
  ! loop over Gauss-points 
   
  do ii=1,numgausspoints
   
    ! calculate shape-function values with respect to natural coordinates        ! == delta3

    call quad8_ansatzfunction(xi_vec(ii), eta_vec(ii), quad8_node_coords_lo, .true., Ni, dNidx, dNidy, detJac)
    
    ! Calculate shape function derivatives with respect to local cartesian coordinate system
    
    do jj=1,8
       
      ! fill operator-matrices
      
      Bu(1,jj) = dNidx(jj)
      Bu(2,jj) = 0.D0
      Bu(3,jj) = dNidy(jj)
      
      Bv(1,jj) = 0.D0
      Bv(2,jj) = dNidy(jj)
      Bv(3,jj) = dNidx(jj)
      
      if (int_selective == .true.) then
         
         Bw(1,jj) = dNidx(jj)
         Bw(2,jj) = dNidy(jj)
         
      end if
    
    end do
    
    ! Calculate parts of the stiffness matrix with full gauss-quadrature
    
    ! Membrane part of element stiffness matrix
    
    facDP = w_xi(ii)*w_eta(ii)*detJac
    
    fac8t3 = matmul(transpose(Bu),Amatel)
    
    kuu = kuu+matmul(fac8t3,Bu)*facDP
    kuv = kuv+matmul(fac8t3,Bv)*facDP
    kvv = kvv+matmul(matmul(transpose(Bv),Amatel),Bv)*facDP
    
    ! Membrane-bending-coupling part of element stiffness matrix
    
    fac8t3 = matmul(transpose(Bu),Bmatel)
    
    kubex = kubex-matmul(fac8t3,Bv)*facDP
    kubey = kubey+matmul(fac8t3,Bu)*facDP
    kvbex = kvbex-matmul(matmul(transpose(Bv),Bmatel),Bv)*facDP
    kvbey = kvbey+matmul(matmul(transpose(Bv),Bmatel),Bu)*facDP
    
    ! Bending part of element stiffness matrix
    
    fac8t3 = matmul(transpose(Bv),Dmatel)
    
    kbexbex = kbexbex+matmul(fac8t3,Bv)*facDP
    kbexbey = kbexbey-matmul(fac8t3,Bu)*facDP
    kbeybey = kbeybey+matmul(matmul(transpose(Bu),Dmatel),Bu)*facDP
    
    ! calculate bending part of stiffness matrix with selective reduced integration
    
    if (int_selective == .true.) then
       
      kww = kww+matmul(matmul(transpose(Bw),Hmatel),Bw)*facDP
       
    end if
    
    ! Penalty stiffness - approach according to Zienkiewicz
    
    SOmegaT = SOmega*k6rot
    
    facDP = 0.25D0*facDP*SOmegaT
    
    do jj=1,8
      do kk=1,8
        ! Diese Werte kommen aus der Penalty-Steifigkeit. Sie ist gedacht, um für Rotz-Rotz eine sinnvolle Steifigkeit zu haben.
        ! In dieser Formulieren werden aber auch alle anderen Steifigkeiten beeinflusst. Dies führt zu einer schlecht konditionierten Matrix für große Probleme.
        ! Aus diesem Grund werden diese Terme weggelassen. Dies wird auch in Ansys bei Shell99 so gehandhabt und ist überrpüft.
        !kuu(jj,kk) = kuu(jj,kk)+dNidy(jj)*dNidy(kk)*facDP
        !kuv(jj,kk) = kuv(jj,kk)-dNidy(jj)*dNidx(kk)*facDP
        !kvv(jj,kk) = kvv(jj,kk)+dNidx(jj)*dNidx(kk)*facDP
        
        !kubez(jj,kk) = kubez(jj,kk)+2.D0*dNidy(jj)*Ni(kk)*facDP
        !kvbez(jj,kk) = kvbez(jj,kk)-2.D0*dNidx(jj)*Ni(kk)*facDP
        
        kbezbez(jj,kk) = kbezbez(jj,kk)+4.D0*Ni(jj)*Ni(kk)*facDP
        
      end do
    end do
   
  end do
  
  area = detJac*4.d0
  
! ******* Anpassung um Ansys-Ergebnisse zu erhalten *******************
  facDP = 1.2d0/max(1.2d0, 1.d0+0.2d0*area/(25.d0*thick**2))
  if (int_selective == .true.) then
    kww = kww*facDP
  end if
  Hmatelt = Hmatel*facDP
! *********************************************************************
  
!  call write_quad_matrix(kbezbez,8,'kbezbez = ')

  deallocate(xi_vec, eta_vec, w_xi, w_eta)

  ! Transverse shear part of element stiffness matrix

  numgausspoints = inplaneintppdts*inplaneintppdts                                                    ! total number of Gauss Points
  
  allocate(xi_vec(numgausspoints), eta_vec(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))
  
  call gauss_integration_values (inplaneintppdts, xi_vec, eta_vec, w_xi, w_eta)
   
  do ii=1,numgausspoints
   
    ! calculate shape-function values with respect to natural coordinates        ! == delta3
    
    call quad8_ansatzfunction(xi_vec(ii), eta_vec(ii), quad8_node_coords_lo, .true., Ni, dNidx, dNidy, detJac)
    
    ! Calculate shape function derivatives with respect to local cartesian coordinate system
    
    do jj=1,8
       
       ! fill operator-matrices
       
       Bbex(1,jj) = 0.D0
       Bbex(2,jj) = -Ni(jj)
       
       Bbey(1,jj) = Ni(jj)
       Bbey(2,jj) = 0.D0
                
       Bw(1,jj)   = dNidx(jj)
       Bw(2,jj)   = dNidy(jj)
       
    end do
    
    facDP=w_xi(ii)*w_eta(ii)*detJac
    
    !write(*,*) 'detJac = ', detJac
    
    ! if reduced integration for bending part of element stiffness
    
    if (int_selective == .false.) then
       
       kww = kww+matmul(matmul(transpose(Bw),Hmatelt),Bw)*facDP
       
    end if
    
    fac8t2 = matmul(transpose(Bw),Hmatelt)
    
    kwbex = kwbex+matmul(fac8t2,Bbex)*facDP
    kwbey = kwbey+matmul(fac8t2,Bbey)*facDP
    
    fac8t2 = matmul(transpose(Bbex),Hmatelt)
    
    kbexbex = kbexbex+matmul(fac8t2,Bbex)*facDP
    kbexbey = kbexbey+matmul(fac8t2,Bbey)*facDP
    
    kbeybey = kbeybey+matmul(matmul(transpose(Bbey),Hmatelt),Bbey)*facDP
   
  end do

  deallocate(xi_vec, eta_vec, w_xi, w_eta)
  
  ! Wenn die Originalwerte belassen werden, ist die Matrix schlecht konditioniert.
  ! Nach dem Vorgehen in ANSYS wird der (1,1)-Wert von kbezbez genommen und auf der
  ! Hauptdiagonalen gesetzt. Die Terme der oberen Dreiecksmatrix entsprechen in ANSYS
  ! minus ein siebtel dieses Wertes. Dieses Vorgehen wird hier kopiert.
  ! Allerdings stimmen die Werte nicht überein.
  do jj=1,8
    do kk=1,8
 
      if (jj .eq. kk) then
        kbezbez(jj,kk) = kbezbez(1,1)
      else
        kbezbez(jj,kk) = -kbezbez(1,1)/7.d0
      end if  
   
    end do
  end do

  ! sort for nodal degrees of freedom
  ! (u1,u2,u3,u4,u5,u6,u7,u8,v1,v2,v3,v4,...,bez8) -> (u1,v1,w1,bex1,bey1,bez1,...,bez8)  

  call sort_shell_dof (8,48,kuu,kuv,kubex,kubey,kubez,kvv,kvbex,kvbey,kvbez,kww,kwbex,kwbey,kbexbex,kbexbey,kbeybey,kbezbez,quad8_K)

  ! Transformation local to global coordinates

  quad8_K = matmul(matmul(TRMat,quad8_K),transpose(TRMat))

!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)  		    'quad8_stiff_mod'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

  return

end subroutine quad8_stiff
