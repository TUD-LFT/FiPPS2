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
!> control subroutine for obtaining the element stiffness matrix in global
!> coordinates for a 2-node 3D beam element with capabilities for tension/
!> compression, bending and torsion
!
!> @details
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 22.06.2010
!
! =================================================================================================
subroutine beam2_stiff_control (fesim, KaaS)

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
type(fe_simulation)                     :: fesim   !< FE-simulation data
!
! Output
!

!
! Input + Output
!
Mat                                     :: KaaS    !< updated aa-part of global stiffness matrix
!
! Inner
!
double precision, dimension(2,3)        :: beam2_node_coords  ! nodal coordinates of element nodes in global coordinates
integer, dimension(3)                   :: beam2_node_ids
double precision, dimension(12,12)      :: beam2_K            ! beam2-element stiffness matrix
double precision                        :: AA                 ! beam cross-section area
double precision                        :: I11                ! beam cross-section moment of inertia about z = Izz
double precision                        :: I22                ! beam cross-section moment of inertia about y = Iyy
double precision                        :: I12                ! beam cross-section product of inertia = Izy
double precision                        :: It                 ! beam cross-section torsional moment of inertia
double precision                        :: Emat               ! material youngs modulus
double precision                        :: Gmat               ! material shear modulus
double precision                        :: lengthWeight, length

integer                                 :: ii,jj
integer                                 :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
beam2_node_coords=0.D0
beam2_K=0.D0
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
      
   ! get beam2-element material stiffness and section properties
   
   call beam2_get_properties(fesim,ii,AA,I11,I22,I12,It,Emat,Gmat)
   
   ! calculate element stiffnes matrix in global coordinates
   
   call beam2_stiff(fesim,ii,beam2_node_coords,AA,I11,I22,It,Emat,Gmat,beam2_K,length)

   ! assemble element stiffness matrix in total stiffness matrix
 
   !call assemble(2,beam2_node_ids,beam2_K,KaaS,KabS,KbbS,dim_dof,num_dof_vec)
   call assemble_aa(fesim,2,beam2_node_ids,beam2_K,KaaS)

   lengthWeight = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%lengthWeight
   
   fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%weight = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%weight + length*lengthWeight

end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_stiff_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine beam2_stiff_control
