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
!> control subroutine for obtaining the element geometric stiffness matrix in
!> global coordinates for a 2-node 3D beam element
!
!> @details
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 24.06.2010
!
!> $Id: beam2_geostiff_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_geostiff_control (fesim,dispg,ntemps,etemps,KgaaS)
!
use konstanten
use fesimulation_typen
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
type(fe_simulation)                             :: fesim                !< FE-Simulation data
double precision, dimension(:), intent(in)      :: dispg(fesim%num_dof) !< degrees of freedom vector with applied boundary conditions
double precision, dimension(fesim%num_nodes), intent(in)      :: ntemps !< (num_dof)-Array containing the nodal temperatures
double precision, dimension(size(fesim%elemente%beam2s,1)), intent(in)   :: etemps !< (num_elements)-Array containing the elemental temperatures
!
! Output
!

!
! Input + Output
!
Mat                                             :: KgaaS                !< updated aa-part of global geometric stiffness matrix
!
! Inner
!
double precision, dimension(2,3)                :: beam2_node_coords    ! nodal coordinates of element nodes in global coordinates
integer, dimension(2)                           :: beam2_node_ids
double precision, dimension(12)                 :: disp_beam2
double precision, dimension(12,12)              :: beam2_Kg             ! beam2-element stiffness matrix
double precision                                :: AA                   ! beam cross-section area
double precision                                :: I11                  ! beam cross-section moment of inertia about
double precision                                :: I22                  ! beam cross-section moment of inertia about
double precision                                :: I12                  ! beam cross-section product of inertia
double precision                                :: It                   ! beam cross-section torsional moment of inertia
double precision                                :: Emat                 ! material youngs modulus
double precision                                :: Gmat                 ! material shear modulus
double precision                                :: btemp                ! beam temperature
double precision                                :: alpha                ! thermal expansion coefficient

integer                                         :: pos,ii,jj,kk
integer                                         :: err_code=0

! character(20)                                   :: name
!
! =================================================================================================
!
! Initialisation
!
beam2_node_coords=0.D0
beam2_Kg=0.D0
!
! =================================================================================================
!
! Calculation
!
! Loop over all 2-node 3D beam elements
do ii=1, size(fesim%elemente%beam2s,1)
   
   ! get node locations
   
   do jj=1,2
      
      beam2_node_ids(jj) = fesim%elemente%beam2s(ii)%nids(jj)
      beam2_node_coords(jj,1:3) = fesim%knoten%nodes(fesim%elemente%beam2s(ii)%nids(jj))%coords(1:3)
      
   end do
   
   ! get global displacement vector for element nodes

   do jj=1,2
      
      pos = (fesim%elemente%beam2s(ii)%nids(jj)-1)*ndof+1
      
      do kk=1,ndof
        disp_beam2(kk+(jj-1)*ndof) = dispg(pos+kk-1)
      end do
      
   end do
    
   ! Get Temperature at element nodes
   ! Set zero, if no temperature is prescribed
   if (fesim%is_beam2temp == .true.) then
     btemp = etemps(ii)
   else if (fesim%is_temperature == .true.) then
     btemp = 0.d0
     do jj = 1,2
       btemp = btemp + ntemps(fesim%elemente%beam2s(ii)%nids(jj))
     end do
     btemp = btemp/2.d0
   end if

   ! get beam2-element thermal expansion coefficient
   alpha = fesim%materialien%mat1s(fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%intMat1ID)%ath
      
   ! get beam2-element material stiffness and section properties
   
   call beam2_get_properties(fesim,ii,AA,I11,I22,i12,It,Emat,Gmat)
   
   ! calculate element geometric stiffnes matrix in global coordinates
   
   call beam2_geostiff(fesim,ii,disp_beam2,beam2_node_coords,AA,I11,I22,Emat,btemp,alpha,beam2_Kg)

   ! assemble element geometric stiffness matrix in total stiffness matrix
   
   call assemble_aa (fesim,2,beam2_node_ids,beam2_Kg, KgaaS)

!    ! write matrix to standard out
!    write(name,'(A4,I16)') 'Kgel',jj
! !    call write_quad_matrix(lsolid20_Kg, 60, name, 'FiPPS')
!    call write_quad_matrix(beam2_Kg, 12, 'Kgel                ', 'ANSYS')
   
end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_geostiff_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine beam2_geostiff_control
