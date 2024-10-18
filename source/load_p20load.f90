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
!> Calculates the equivalent load vector for surface pressure on the 20-node layered
!> solid element Lsolid20
!
!> @details
!> Calculation of the equivalent load vector for the surface pressure on the 20-node
!> layered solid element Lsolid20; the pressure can be applied on any of the 6
!> element surfaces and acts normal on the element surface at the specific point;
!> positive pressure acts towards the inside of the element; the equivalent load
!> vector is added to the global load vector at the corresponding node entries;
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: load_p20load.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine load_p20load (fesim,hh,lvec,pos,pElem)
! =================================================================================================
! Use
!
use mat_func
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
! Input
!
type(fe_simulation)                                           :: fesim
integer, intent(in)                                           :: hh !< Index of current load
integer, intent(in)                                           :: pos !< Position of loadvector in loadvector array
!
! Input + Output
!
double precision, dimension(:,:), intent(inout) :: lvec(1:fesim%num_dof,1:fesim%num_loadcases) !< Total loadvector for current load
double precision, dimension(fesim%num_elements),intent(inout) :: pElem !< Vector of Element pressures for vtk-output
!
! Internal
!
integer                                                       :: loadI, eid
integer                                                       :: err_code=0
!
! =================================================================================================
!
! Calculation
!
 ! loop over p20loads
 do loadI=1,size(fesim%lasten%p20loads,1)

   if (fesim%lasten%loads(hh)%lidi == fesim%lasten%p20loads(loadI)%lid) then

    ! if p20load on one element only
    
    if (fesim%lasten%p20loads(loadI)%thru == .false.) then
       
       ! apply pressure
       call apply_nodal_forces_solid20(fesim%lasten%p20loads(loadI)%eid1, loadI, lvec, pos, fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac)
       
       pElem(fesim%elemente%lsolid20s(fesim%lasten%p20loads(loadI)%eid1)%eid) = pElem(fesim%elemente%lsolid20s(fesim%lasten%p20loads(loadI)%eid1)%eid) - fesim%lasten%p20loads(loadI)%p*fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
    
    ! pressure on multiple elements
    
    else if (fesim%lasten%p20loads(loadI)%thru == .true.) then
      
      ! loop over elements to which pressure shall be applied
      do eid = fesim%lasten%p20loads(loadI)%eid1,fesim%lasten%p20loads(loadI)%eid2

        ! apply pressure
        call apply_nodal_forces_solid20(eid, loadI, lvec, pos, fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac)
        
        pElem(fesim%elemente%lsolid20s(eid)%eid) =  pElem(fesim%elemente%lsolid20s(eid)%eid) - fesim%lasten%p20loads(loadI)%p*fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
      
      end do
       
    else
       
       write(*,*) 'Error processing p20load ',hh
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
   write(*,*)                      'load_p20load'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

contains

  subroutine apply_nodal_forces_solid20 (eid, loadID, lvec, pos, loadFac)
  
  use integration_schemes
  
  implicit none
  !
  ! =================================================================================================
  !
  ! Interface
  !
  interface
    subroutine hexa20_shapefunctions_xi_eta_zeta(xi, eta, zeta, Ni, dNidxi, dNideta, dNidzeta)
      double precision, intent(in)                  :: xi, eta, zeta
      double precision, dimension(20), optional, intent(out)        :: Ni
      double precision, dimension(20), optional, intent(out)        :: dNidxi, dNideta, dNidzeta
    end subroutine hexa20_shapefunctions_xi_eta_zeta
  end interface 
  !
  ! =================================================================================================
  ! 
  ! Data types
  !
  ! Input
  integer, intent(in)                             :: loadID, eid
  integer, intent(in)                             :: pos
  double precision, intent(in)                    :: loadFac
  !
  ! Input + Output
  !
  double precision, dimension(:,:), intent(inout) :: lvec(1:fesim%num_dof,1:fesim%num_loadcases)
  !
  ! Internal
  !
  integer                                         :: surf
  integer                                         :: jj, kk
  integer                                         :: nip, N, ind_begin, ind_end
  integer                                         :: lvec_begin, lvec_end
  integer                                         :: nodeID
  double precision, dimension(:), allocatable     :: xi, eta, zeta
  double precision, dimension(3)                  :: u1, u2, u3, g
  double precision                                :: dir, p_val, weight
  double precision, dimension(3,3)                :: transMat
  double precision, dimension(:), allocatable     :: ip_vec_1, ip_vec_2, w_1, w_2
  double precision, dimension(60)                 :: f_vec, NTg
  double precision, dimension(20,3)               :: node_coords
  double precision, dimension(20)                 :: Ni, dNidxi, dNideta, dNidzeta
  !
  ! =================================================================================================
  !
  ! Initialisation
  !
  nip   = 3     ! Number of integration points in plane used for integration
                ! Standard value for compliance with ANSYS results: nip = 3
  f_vec = 0.d0
  !
  ! =================================================================================================
  !
  ! Calculation
  ! 
  ! Number of integration points
  N = nip*nip
  
  ! allocate arrays
  allocate(xi(N),eta(N),zeta(N),ip_vec_1(N),ip_vec_2(N),w_1(N),w_2(N))

  ! get node coordinates  
  do jj = 1,20
    node_coords(jj,:) = fesim%knoten%nodes(fesim%elemente%lsolid20s(eid)%nids(jj))%coords(1:3)
  end do
  
  ! get ID of surface
  surf = fesim%lasten%p20loads(loadID)%surf
  
  ! get pressure value
  p_val  = fesim%lasten%p20loads(loadID)%p * loadFac
  
  ! get integration point coordinates and weights
  call integration_points_2d(nip, "Gauss", ip_vec_1, w_1, &
                             nip, "Gauss", ip_vec_2, w_2)

  ! set up integration point coordinates according to surface
  ! dir =  1.d0 if surface normal points in same direction as the 
  !             corresponding natural coordinate
  ! dir = -1.d0 if surface normal points in opposite direction as the
  !             corresponding natural coordinate
  if (surf .EQ. 1) then
    xi   = ip_vec_1
    eta  = ip_vec_2
    zeta = -1.d0
    dir  = -1.d0
  else if (surf .EQ. 2) then
    xi   = ip_vec_1
    eta  = -1.d0
    zeta = ip_vec_2
    dir  = -1.d0
  else if (surf .EQ. 3) then
    xi   = 1.d0
    eta  = ip_vec_1
    zeta = ip_vec_2
    dir  = 1.d0
  else if (surf .EQ. 4) then
    xi   = ip_vec_1
    eta  = 1.d0
    zeta = ip_vec_2
    dir  = 1.d0
  else if (surf .EQ. 5) then
    xi   = -1.d0
    eta  = ip_vec_1
    zeta = ip_vec_2
    dir  = -1.d0
  else if (surf .EQ. 6) then
    xi   = ip_vec_1
    eta  = ip_vec_2
    zeta = 1.d0
    dir  = 1.d0
  else
    write(*,*) 'Error processing p20load ',hh
    write(*,*) 'Invalid surface number, valid numbers are:'
    write(*,*) '  1, 2, 3, 4, 5 or 6'
    err_code = 1
    goto 9999
  end if
  
  ! loop over integration points
  do jj = 1,N

    ! numerical integration weight
    weight = w_1(jj) * w_2(jj)

    ! get values of shapefunctions and their derivates
    call hexa20_shapefunctions_xi_eta_zeta(xi(jj),eta(jj),zeta(jj),Ni=Ni,dNidxi=dNidxi,dNideta=dNideta,dNidzeta=dNidzeta)

    ! initialise arrays
    u1   = 0.d0
    u2   = 0.d0
    u3   = 0.d0

    ! get vectors u1, u2, u3 aligned along natural coordinates xi, eta, zeta
    do kk = 1,20
      u1(1) = u1(1) + dNidxi(kk)   * node_coords(kk,1)
      u1(2) = u1(2) + dNidxi(kk)   * node_coords(kk,2)
      u1(3) = u1(3) + dNidxi(kk)   * node_coords(kk,3)
      u2(1) = u2(1) + dNideta(kk)  * node_coords(kk,1)
      u2(2) = u2(2) + dNideta(kk)  * node_coords(kk,2)
      u2(3) = u2(3) + dNideta(kk)  * node_coords(kk,3)
      u3(1) = u3(1) + dNidzeta(kk) * node_coords(kk,1)
      u3(2) = u3(2) + dNidzeta(kk) * node_coords(kk,2)
      u3(3) = u3(3) + dNidzeta(kk) * node_coords(kk,3)
    end do

    ! get surface normal vector g according to surface
    if ((surf .EQ. 1) .OR. (surf .EQ. 6)) then
      g(:) = cross_product3(u1(:),u2(:))
    else if ((surf .EQ. 2) .OR. (surf .EQ. 4)) then
      g(:) = cross_product3(u3(:),u1(:))
    else if ((surf .EQ. 3) .OR. (surf .EQ. 5)) then
      g(:) = cross_product3(u2(:),u3(:))
    end if

    ! turn g to opposite direction, if necessary
    g(:) = dir * g(:)

    ! multiplicate shape function values with vector g
    do kk = 1,20
      NTg((kk-1)*3+1) = Ni(kk) * g(1)
      NTg((kk-1)*3+2) = Ni(kk) * g(2)
      NTg((kk-1)*3+3) = Ni(kk) * g(3)
    end do

    ! integrate
    f_vec = f_vec - NTg * p_val * weight
    
  end do
  
  ! deallocate arrays
  deallocate(xi,eta,zeta,ip_vec_1,ip_vec_2,w_1,w_2)
 
  ! add equivalent nodal forces to global force vector
  do jj = 1,20
    nodeID     = fesim%elemente%lsolid20s(eid)%nids(jj)
    
    ind_begin  = (jj-1)*3+1
    ind_end    = ind_begin+2
    
    lvec_begin = (nodeID-1)*ndof+1
    lvec_end   = lvec_begin+2
    
    if (fesim%knoten%nodes(nodeID)%cid .NE. 0) then
      transMat(1:3,1:3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nodeID)%cid)%transMat(1:3,1:3)
      transMat = transpose(transMat)
      f_vec(ind_begin:ind_end) = matmul(transMat,f_vec(ind_begin:ind_end))
    end if
    
    lvec(lvec_begin:lvec_end,pos) = lvec(lvec_begin:lvec_end,pos)+f_vec(ind_begin:ind_end)
  end do  
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'load_p20load'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if

  end subroutine apply_nodal_forces_solid20
 
end subroutine load_p20load
    
