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
!> Subroutine for computing load-vectors entrys resulting from line load on beam2
!
!> @details
!> Subroutine calculates contribution to the loadvector resulting from line load
!> on all beam elements (only beam2 yet);
!
!> @author Florian Dexl, TU Dresden, WiMi 2017
!
!> $Id: load_p2load.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine load_p2load (fesim,hh,lvec,pos,pElem)
!
! =================================================================================================
!
! Use
!
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
type(fe_simulation)                                           :: fesim
double precision, dimension(:,:), intent(inout)               :: lvec(1:fesim%num_dof,1:fesim%num_loadcases) !< Total loadvector for current load
double precision, dimension(fesim%num_elements),intent(inout) :: pElem !< Vector of Element pressures for vtk-output
integer, intent(in)                                           :: pos !< Position of loadvector in loadvector array
integer, intent(in)                                           :: hh !< Index of current load

integer                                                       :: loadI, eid
integer                                                       :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
!
! =================================================================================================
!
! Calculation
!
do loadI=1,size(fesim%lasten%p2loads,1)
  if (fesim%lasten%loads(hh)%lidi == fesim%lasten%p2loads(loadI)%lid) then
     
    ! if p2load on one element only
    
    if (fesim%lasten%p2loads(loadI)%thru == .false.) then
       
       call apply_nodal_forces_beam2(fesim%lasten%p2loads(loadI)%eid1, loadI, fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac)
       
       pElem(fesim%elemente%beam2s(fesim%lasten%p2loads(loadI)%eid1)%eid) = pElem(fesim%elemente%beam2s(fesim%lasten%p2loads(loadI)%eid1)%eid) + fesim%lasten%p2loads(loadI)%pi(1)*fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
    
    ! p2 load on multiple elements
    
    else if (fesim%lasten%p2loads(loadI)%thru == .true.) then
    
      do eid=fesim%lasten%p2loads(loadI)%eid1,fesim%lasten%p2loads(loadI)%eid2
      
        call apply_nodal_forces_beam2(eid, loadI, fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac)
        
        pElem(fesim%elemente%beam2s(eid)%eid) =  pElem(fesim%elemente%beam2s(eid)%eid) + fesim%lasten%p2loads(loadI)%pi(1)*fesim%lasten%loads(hh)%sfaci*fesim%lasten%loads(hh)%sfac
      
      end do
       
    else
       
       write(*,*) 'Error processing p2load ',hh
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
   write(*,*)                      'load_p2load'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

contains

  subroutine apply_nodal_forces_beam2 (eid, loadID, loadFac)
  
  use mat_func
  
  implicit none
  !
  ! =================================================================================================
  !
  ! Interface
  !
  !
  ! =================================================================================================
  ! 
  ! Data types
  !  
  integer, intent(in)                               :: loadID, eid
  double precision, intent(in)                      :: loadFac

  integer                                           :: kk
  integer                                           :: index, nid
  integer, dimension(2)                             :: node_ids_beam2
  
  double precision                                  :: elem_length
  double precision                                  :: force_val, moment_val
  double precision, dimension(2,3)                  :: force_vec, moment_vec
  double precision, dimension(3)                    :: vv
  double precision, dimension(2,3)                  :: node_coords
  double precision, dimension(3,3)                  :: transMat
  double precision, dimension(12,12)                :: TRMat
  !
  ! =================================================================================================
  !
  ! Calculation
  !
         
  ! get beam2-element node-ids
  node_ids_beam2(1:2) = fesim%elemente%beam2s(eid)%nids(1:2)
  
  ! search for node-coordinates
  do kk = 1,2
    node_coords(kk,1:3) = fesim%knoten%nodes(fesim%elemente%beam2s(eid)%nids(kk))%coords(1:3)
  end do
  
  ! get transformation matrix from local element coordinate system
  ! to global coordinate system and element length
  vv(1:3)=fesim%elemente%beam2s(eid)%xi(1:3)
  call beam2_rotation (node_coords,vv,'lg',TRMat,elem_length)
  
  ! Vektoren der Kraefte und Momente in lokalen Elementkoordinaten
  if (fesim%lasten%p2loads(loadID)%dir .EQ. 1) then
    ! Positiver Wert in negative z-Richtung
    force_vec(:,1)  =  0.d0
    force_vec(:,2)  =  0.d0
    force_vec(:,3)  = -1.d0
    moment_vec(:,1) =  0.d0
    moment_vec(:,2) =  1.d0
    moment_vec(:,3) =  0.d0
  else if (fesim%lasten%p2loads(loadID)%dir .EQ. 2) then
    ! Positiver Wert in negative y-Richtung
    force_vec(:,1)  =  0.d0
    force_vec(:,2)  = -1.d0
    force_vec(:,3)  =  0.d0
    moment_vec(:,1) =  0.d0
    moment_vec(:,2) =  0.d0
    moment_vec(:,3) = -1.d0
  end if    
  
  ! Berechnung der aequivalenten Knotenkraft und -moment
  force_val  = loadFac*fesim%lasten%p2loads(loadID)%pi(1)*elem_length/2.d0
  moment_val = loadFac*fesim%lasten%p2loads(loadID)%pi(1)*elem_length**2.d0/12.d0
  force_vec(:,:)  = force_val*force_vec(:,:)
  moment_vec(1,:) = moment_val*moment_vec(1,:)
  moment_vec(2,:) = -1.d0*moment_vec(1,:)
  
  ! Addieren der Kraefte und Momente zu Lastvektor
  do kk = 1, 2                                                              ! Schleife über alle 2 Knoten
      
    ! Transformieren der Kraft- und Momentenvektoren in globales Koordinatensystem
    force_vec(kk,:) = matmul(TRMat(1:3,1:3),force_vec(kk,:))
    moment_vec(kk,:) = matmul(TRMat(1:3,1:3),moment_vec(kk,:))
      
    nid = node_ids_beam2(kk)
    
    if (fesim%knoten%nodes(nid)%cid .ne. 0) then
      transMat(1:3,1:3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(nid)%cid)%transMat(1:3,1:3)
      transMat = transpose(transMat)
      force_vec(kk,:) = matmul(transMat,force_vec(kk,:))
      moment_vec(kk,:) = matmul(transMat,moment_vec(kk,:))
    end if
  
    index = (nid-1)*ndof+1
    
    lvec(index:index+2,pos)   = lvec(index:index+2,pos) + force_vec(kk,:)
    lvec(index+3:index+5,pos) = lvec(index+3:index+5,pos) + moment_vec(kk,:)
    
  end do
  
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'load_p2load'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if

  end subroutine apply_nodal_forces_beam2
 
end subroutine load_p2load
