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
!> subroutine for computing load-vectors entrys resulting from surface
!> pressure
!
!> @details
!> subroutine calculates contribution to the loadvector resulting from surface
!> pressure on all 2D quadrilateral elements (only quad8 yet)
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter, 24.10.2012
!
!> $Id: load_p8load.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine load_p8load (fesim,hh,lvec,pos,pElem)

use konstanten
use fesimulation_typen
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
type(fe_simulation)                             :: fesim
double precision, dimension(:,:), intent(inout) :: lvec(1:fesim%num_dof,1:fesim%num_loadcases)      !< loadvectors for different loads
integer, intent(in)                             :: pos      !< column of current load
integer, intent(in)                             :: hh       !< loop index over loads

integer                                         :: inplaneintppd
integer                                         :: loadI, eid

integer                                         :: err_code=0

double precision, dimension(fesim%num_elements) :: pElem
!
! =================================================================================================
!
! Initialisation
!
inplaneintppd = 2                                ! Anzahl der Integrationspunkte über der Elementmittelfläche
!
! =================================================================================================
!
! Calculation
!

do loadI=1,size(fesim%lasten%p8loads,1)
  if (fesim%lasten%loads(hh)%lidi == fesim%lasten%p8loads(loadI)%lid) then
     
    ! if p8load on one element only
    
    if (fesim%lasten%p8loads(loadI)%thru == .false.) then
       
       call apply_nodal_forces_quad8(fesim%lasten%p8loads(loadI)%eid1, loadI, fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac)
        
       pElem(fesim%elemente%quad8s(fesim%lasten%p8loads(loadI)%eid1)%eid) = pElem(fesim%elemente%quad8s(fesim%lasten%p8loads(loadI)%eid1)%eid) - fesim%lasten%p8loads(loadI)%pi(1)*fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
    
    ! pressure on multiple elements
    
    else if (fesim%lasten%p8loads(loadI)%thru == .true.) then
    
      do eid=fesim%lasten%p8loads(loadI)%eid1,fesim%lasten%p8loads(loadI)%eid2
      
        call apply_nodal_forces_quad8(eid, loadI, fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac)
        
        pElem(fesim%elemente%quad8s(eid)%eid) =  pElem(fesim%elemente%quad8s(eid)%eid) - fesim%lasten%p8loads(loadI)%pi(1)*fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
      
      end do
       
    else
       
       write(*,*) 'Error processing p8load ',hh
       err_code = 1
       goto 9998
       
    end if
  end if
end do

!
! =================================================================================================
!
! Error handling
!
9998 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'load_p8load'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

contains

  subroutine apply_nodal_forces_quad8 (eid, loadID, loadFac)
  
  use integration_schemes
  
  implicit none
  !
  ! =================================================================================================
  !
  ! Interface
  !
    interface 
      subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
        double precision, intent(in)                                :: xi, eta
        
        double precision, dimension(8), optional, intent(out)       :: Ni
        double precision, dimension(8), optional, intent(out)       :: dNidxi, dNideta
        double precision, dimension(8), optional, intent(out)       :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
      end subroutine quad8_ansatzfunction_xieta
    end interface
  !
  ! =================================================================================================
  ! 
  ! Data types
  !
  integer, intent(in)                               :: loadID, eid
  double precision, intent(in)                      :: loadFac
  
  integer                                           :: numgausspoints
  integer                                           :: kk, mm
  integer                                           :: index, nid
  integer, dimension(8)                             :: node_ids_quad8
  
  double precision, dimension(3)                    :: force
  double precision, dimension(3)                    :: nor_vec
  double precision, dimension(:), allocatable       :: xi_int, eta_int, w_xi, w_eta
  double precision, dimension(8,3)                  :: node_coords
  double precision, dimension(8)                    :: Ni, dNidxi, dNideta
  double precision                                  :: dxdxi, dxdeta, dydxi, dydeta, dzdxi, dzdeta
  double precision                                  :: Exieta, Fxieta, Gxieta
  double precision, dimension(3,3)                  :: transMat
  !
  ! =================================================================================================
  !
  ! Calculation
  !
         
  ! get quad8-element node-ids
  
  node_ids_quad8(1:8) = fesim%elemente%quad8s(eid)%nids(1:8)
  
  ! search for node-coordinates
  
  do kk = 1,8
    node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%quad8s(eid)%nids(kk))%coords(1:3)
  end do
  
  numgausspoints = inplaneintppd*inplaneintppd                              ! total number of Gauss Points
  
  allocate(xi_int(numgausspoints), eta_int(numgausspoints))
  allocate(  w_xi(numgausspoints),   w_eta(numgausspoints))
  
  call integration_points_2d(inplaneintppd, "Gauss",  xi_int, w_xi, &
                             inplaneintppd, "Gauss", eta_int, w_eta)
  
  do kk = 1, 8                                                              ! Schleife über alle 8 Knoten und deren Ansatzfunktionen
  
    force = 0.d0
  
    ! Integration der Ansatzfunktion des kk-ten Knotens und Berechung der arbeitsäquivalenten Kraft
  
    do mm = 1, numgausspoints                                               ! Integration der Ansatzfunktion den entsprechenden Knotens
    
      call quad8_ansatzfunction_xieta(xi_int(mm), eta_int(mm), Ni=Ni, dNidxi=dNidxi, dNideta=dNideta)
  
      ! Calculate Jacobian
      
      dxdxi  = SUM(dNidxi(1:8) *node_coords(1:8,1))
      dxdeta = SUM(dNideta(1:8)*node_coords(1:8,1))
      
      dydxi  = SUM(dNidxi(1:8) *node_coords(1:8,2))
      dydeta = SUM(dNideta(1:8)*node_coords(1:8,2))
      
      dzdxi  = SUM(dNidxi(1:8) *node_coords(1:8,3))
      dzdeta = SUM(dNideta(1:8)*node_coords(1:8,3))
      
      Exieta = dxdxi**2 + dydxi**2 + dzdxi**2
      Fxieta = dxdxi*dxdeta + dydxi*dydeta + dzdxi*dzdeta
      Gxieta = dxdeta**2 + dydeta**2 + dzdeta**2
    
      if (fesim%lasten%p8loads(loadI)%cid == 0) then
        ! Berechnung der Normalen auf die Elementmittelsfläche am Integrationspunkt
        nor_vec(1) = dydxi*dzdeta - dzdxi*dydeta
        nor_vec(2) = dzdxi*dxdeta - dxdxi*dzdeta
        nor_vec(3) = dxdxi*dydeta - dydxi*dxdeta
    
        nor_vec = nor_vec / SQRT(DOT_PRODUCT(nor_vec,nor_vec))
      else if (fesim%lasten%p8loads(loadI)%cid == -1) then
        ! Orientierung Kraftvektor entlang x-Achse aus globalem Koordinatensystem
        nor_vec(1) = 1.d0
        nor_vec(2) = 0.d0
        nor_vec(3) = 0.d0
      else if (fesim%lasten%p8loads(loadI)%cid == -2) then
        ! Orientierung Kraftvektor entlang y-Achse aus globalem Koordinatensystem
        nor_vec(1) = 0.d0
        nor_vec(2) = 1.d0
        nor_vec(3) = 0.d0
      else if (fesim%lasten%p8loads(loadI)%cid == -3) then
        ! Orientierung Kraftvektor entlang z-Achse aus globalem Koordinatensystem
        nor_vec(1) = 0.d0
        nor_vec(2) = 0.d0
        nor_vec(3) = 1.d0
      else
        ! Orientierung Kraftvektor entlang z-Achse aus nutzerdefiniertem Koordinatensystem
        nor_vec(1:3) = fesim%koordinatensysteme%coords(fesim%lasten%p8loads(loadI)%cid)%transMat(1:3,3)
      end if
  
      force = force + nor_vec*Ni(kk)*w_xi(mm)*w_eta(mm)*SQRT(Exieta*Gxieta-Fxieta**2)
  
    end do
    
    force = force * fesim%lasten%p8loads(loadID)%pi(1) * loadFac
      
    nid = node_ids_quad8(kk)
    
    if (fesim%knoten%nodes(nid)%cid .ne. 0) then
      transMat(1:3,1:3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nid)%cid)%transMat(1:3,1:3)
      transMat = transpose(transMat)
      force = matmul(transMat,force)
    end if
  
    index = (nid-1)*ndof+1
    
    lvec(index:index+2,pos) = lvec(index:index+2,pos) + force
    
  end do
  
  deallocate(xi_int, eta_int)
  deallocate(  w_xi,   w_eta)
  
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'load_p8load'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if

  end subroutine apply_nodal_forces_quad8
 
end subroutine load_p8load
    
