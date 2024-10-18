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
!
!> @details
!
!> @author 
!
!> $Id: sol2_get_geometric_stiffness_matrix.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol2_get_geometric_stiffness_matrix(fesim, scloop, KgaaS, Utot)

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
  type(fe_simulation)                                   :: fesim
  integer,intent(in)                                    :: scloop
  double precision, dimension(fesim%num_dof), intent(in):: Utot
  double precision, allocatable, dimension(:)           :: ntemps, beam2temps, quad8temps, lsolid20temps
  
  integer                                               :: err_code=0
  
  Mat                                                   :: KgaaS
  
  PetscMPIInt                                           :: rank
  PetscErrorCode                                        :: ierr
  integer                                               :: start, ende, ctr, ctm
  
!
! =================================================================================================
!
! Initialisation
!
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
!
! =================================================================================================
!
! Calculation
!  
   allocate(ntemps(fesim%num_nodes), beam2temps(size(fesim%elemente%beam2s,1)), quad8temps(size(fesim%elemente%quad8s,1)), lsolid20temps(size(fesim%elemente%lsolid20s,1)))
   ntemps = 0.d0
   beam2temps = 0.d0
   quad8temps = 0.d0
   lsolid20temps = 0.d0
   ! Get nodal temperatures
   if (fesim%is_temperature == .true.) call get_node_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, ntemps)
   ! Get elemental temperatures at beam2 elements
   if (fesim%is_beam2temp .eqv. .true.) call get_beam2_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, beam2temps)
   ! Get elemental temperatures at quad8 elements
   if (fesim%is_quad8temp .eqv. .true.) call get_quad8_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, quad8temps)
   ! Get elemental temperatures at lsolid20 elements
   if (fesim%is_lsolid20temp .eqv. .true.) call get_lsolid20_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, lsolid20temps)

  if (fesim%is_quad8 == .true.) then
    if (rank == 0) call SYSTEM_CLOCK(start,ctr,ctm)
    call quad8_geostiff_control (fesim, Utot, ntemps, quad8temps, KgaaS)
    if (rank == 0) call SYSTEM_CLOCK(ende,ctr,ctm)
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Quad8-Zeit', (ende-start)/DBLE(ctr)
  end if

  if (fesim%is_lsolid20 == .true.) then
    if (rank == 0) call SYSTEM_CLOCK(start,ctr,ctm)
    call lsolid20_geostiff_control (fesim, scloop, Utot, ntemps, lsolid20temps, KgaaS)
    if (rank == 0) call SYSTEM_CLOCK(ende,ctr,ctm)
    if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) 'Lsolid20-Zeit', (ende-start)/DBLE(ctr)
  end if
  
  if (fesim%is_beam2 == .true.) then
    call beam2_geostiff_control (fesim, Utot, ntemps, beam2temps, KgaaS)
  end if

  deallocate(ntemps, beam2temps, quad8temps, lsolid20temps)
!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol2_get_geometric_stiffness_matrix'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine sol2_get_geometric_stiffness_matrix
