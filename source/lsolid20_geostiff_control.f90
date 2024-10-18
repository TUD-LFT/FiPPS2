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
!> control subroutine for computing element geoemtric stiffness matrix for quadratic 
!> layered 20-node solid element and assemble it to the total
!> stiffness matrix
!
!> @details
!> subroutine calculates the element geometric stiffness matrix for quadratic layered 20-node
!> solid element and assembles it to the total stiffness matrix
!
!> @author Florian Dexl, TU Dresden, wiss. Mitarbeiter, 27.09.2018
!
!> $Id: lsolid20_geostiff_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine lsolid20_geostiff_control (fesim, scloop, dispg, ntemps, etemps, KgaaS)

! =================================================================================================
! Use
!
use konstanten
use fesimulation_typen
use globale_variablen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscmat.h"
use petscmat
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
type(fe_simulation), intent(in)                               :: fesim  !< FE-Simulation data
integer, intent(in)                                           :: scloop !< Number of current subcase
double precision, dimension(fesim%num_dof), intent(in)        :: dispg  !< (num_dof)-Array containing the displacements in the nodal coordinate system
double precision, dimension(fesim%num_nodes), intent(in)      :: ntemps !< (num_dof)-Array containing the nodal temperatures
double precision, dimension(size(fesim%elemente%lsolid20s,1)), intent(in)   :: etemps !< (num_dof)-Array containing the elemental temperatures
!
! Output
!

!
! Input + Output
!
Mat                                             :: KgaaS                    !< AA-Part of of the global geometric stiffness matrix (sparse - PETSC)
!
! Internal
!
integer, dimension(20)              :: node_ids
double precision, dimension(20)     :: nodetemps
double precision, dimension(60)     :: disp_l20
double precision, dimension(20,3)   :: node_coords
double precision, dimension(60,60)  :: lsolid20_Kg
double precision, dimension(120,120):: Kg
double precision, allocatable       :: Clay(:,:,:), ath_lay(:,:), lth(:), initstress_glob(:,:,:)
integer, allocatable, dimension(:)  :: nop

double precision                    :: tth

integer                             :: ii, jj, kk, lay, ind, facI, nlay, nip, cid
integer                             :: err_code=0

! character(20)                       :: name
!
! =================================================================================================
!
! Initialisation
!
 facI   = 1
 nip = 2    ! number of integration points in plane
!
! =================================================================================================
!
! Calculation
!

! Loop over all layered solid elements
 do jj = 1,size(fesim%elemente%lsolid20s,1)

   if (textoutput == .true. .and. jj == facI*100) then
    
     write(*,*) 'aktuelles Lsolid20-Element: ', facI*100
    
     facI = facI+1
    
   end if

   ! get number of layers for current element
   nlay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%lay
   
   ! allocate arrays dependent on layer number
   allocate(Clay(6,6,nlay))
   allocate(ath_lay(6,nlay))
   allocate(lth(nlay))
   allocate(nop(nlay))
   
   ! get stiffness matrices for layers
   Clay = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%C(:,:,:)
   ! get thermal expansion coefficients for layers
   ath_lay(:,:) = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%ath(:,:)
   ! get thickness for layers
   lth  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%lth(:)
   ! get total thickness of laminat
   tth  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%tth
   ! get number of integration points over thickness for each layer
   nop  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%nop(:)
   ! get ID of coordinate system for orientation of element coosy
   cid  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(jj)%int_pid)%cid

   do kk = 1,20
     node_ids(kk)            = fesim%elemente%lsolid20s(jj)%nids(kk)
     ind = (node_ids(kk)-1)*ndof
     do ii = 1,3
       disp_l20(ii+(kk-1)*3) = dispg(ind+ii)
     end do
     node_coords(kk,1:3)     = fesim%knoten%nodes(fesim%elemente%lsolid20s(jj)%nids(kk))%coords(1:3)
   end do

   ! Get temperature at element nodes
   nodetemps = 0.d0
   if (fesim%is_lsolid20temp == .true.) then
     nodetemps(:)    = etemps(jj)
   else if (fesim%is_temperature == .true.) then
     do kk = 1,20
       nodetemps(kk) = ntemps(fesim%elemente%lsolid20s(jj)%nids(kk))
     end do
   end if

   ! Get initial stress in global coordinate system
   allocate(initstress_glob(nlay,6,nip*nip*maxval(nop)))
   initstress_glob(:,:,:) = 0.d0
   if (fesim%is_multistep .EQ. .TRUE.) then
     if (fesim%lasten%subcases(scloop)%upstress .EQ. .TRUE.) then ! if initial stress shall be applied
       do lay = 1,nlay
         initstress_glob(lay,:,1:nip*nip*nop(lay)) = fesim%residuals%initstr20s(jj)%lay(lay)%is_ip(:,:)
       end do
     end if
   end if

   ! calculate element stiffness matrix  
   call lsolid20_geostiff(fesim,scloop,nlay,Clay,ath_lay,lth,tth,node_coords,nip,nop,cid,lsolid20_Kg,disp_l20,nodetemps,initstress_glob)

   ! Blow up stiffness matrix to dof = 6
   ! Zeros are added at entries of rotational dofs
   call blow_up_stiff(20, lsolid20_Kg, Kg)

   call assemble_aa (fesim,20,node_ids(:),Kg,KgaaS)

   deallocate(Clay,ath_lay,lth,nop,initstress_glob)

!    ! write matrix to standard out
!     write(name,'(A4,I16)') 'Kgel',jj
! !     call write_quad_matrix(lsolid20_Kg, 60, name, 'FiPPS')
!     call write_quad_matrix(lsolid20_Kg, 60, 'Kgel                ', 'ANSYS')

 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'lsolid20_geostiff_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine lsolid20_geostiff_control
