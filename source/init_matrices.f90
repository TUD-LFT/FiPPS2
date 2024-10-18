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
!> $Id: init_matrices.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine init_matrices(fesim, KaaS, KgaaS, rank)

  use konstanten
  use fesimulation_typen
  use pre_assemble_types
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

#if defined (PETSC_HAVE_MPIUNI)
    MPIUNI_FInt MPI_INTEGER,MPI_LOGICAL
#endif
!
! =================================================================================================
!
! Interfaces
!
  INTERFACE 
    SUBROUTINE CLEAR_NODE_NODE_CON(VIRTNODES,NUM_NODES)
      USE PRE_ASSEMBLE_TYPES
      INTEGER, INTENT(IN) :: NUM_NODES
      TYPE (VIRTUALNODE) ,TARGET :: VIRTNODES(NUM_NODES)
    END SUBROUTINE CLEAR_NODE_NODE_CON
  END INTERFACE
  INTERFACE 
    SUBROUTINE GET_NODE_NODE_CON(FESIM,VIRTNODES)
      USE FESIMULATION_TYPEN
      USE PRE_ASSEMBLE_TYPES
      TYPE(FE_SIMULATION) :: FESIM
      TYPE (VIRTUALNODE) ,TARGET :: VIRTNODES(FESIM%NUM_NODES)
    END SUBROUTINE GET_NODE_NODE_CON
  END INTERFACE
  INTERFACE 
    SUBROUTINE ESTIMATE_NUMBER_NONZEROS(FESIM,VIRTNODES,RANGES, &
    & MPI_SIZE,D_NZZ,O_NZZ,GEOMETRY)
      USE FESIMULATION_TYPEN
      USE PRE_ASSEMBLE_TYPES
      TYPE(FE_SIMULATION) :: FESIM
      INTEGER(KIND=4), INTENT(IN) :: MPI_SIZE
      TYPE (VIRTUALNODE) ,TARGET, INTENT(IN) :: VIRTNODES(      &
              &FESIM%NUM_NODES)
      INTEGER(KIND=4), INTENT(IN) :: RANGES(MPI_SIZE+1)
      INTEGER(KIND=4), INTENT(OUT) :: D_NZZ(FESIM%INTERNALS%DIM_DOF)
      INTEGER(KIND=4), INTENT(OUT) :: O_NZZ(FESIM%INTERNALS%DIM_DOF)
      LOGICAL(KIND=4), INTENT(IN) :: GEOMETRY
    END SUBROUTINE ESTIMATE_NUMBER_NONZEROS
  END INTERFACE  
!
! =================================================================================================
!
! Data types
!
  type(fe_simulation)   :: fesim
  integer               :: ii, size
  
  integer, dimension(:), allocatable      :: ranges 
  
  type(virtualNode), dimension(:), allocatable, target :: virtNodes
  
  PetscInt       :: m
  PetscInt       :: bs
  PetscInt, dimension(:), allocatable  :: d_nzz
  PetscInt, dimension(:), allocatable  :: o_nzz
  PetscErrorCode :: ierr
  PetscInt                             :: start, ende
  PetscMPIInt, intent(in)              :: rank

  Mat            :: KaaS         ! Steifigkeitsmatrix ohne gesperrte Freiheiten als Sparse-Matrix
  Mat            :: KgaaS
  
  integer       :: err_code=0
    integer                                               :: st, en, ctr, ctm ! Variablen zur Zeitmessung
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
  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  START - Bereite Matrizen vor'
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr); CHKERRQ(ierr)
 
  if (rank .eq. 0) then
    allocate(virtNodes(fesim%num_nodes))
    call get_node_node_con(fesim,virtNodes)
  end if
  
  ! Anlegen von Kaa
  m = fesim%internals%dim_dof
  call MatCreate(PETSC_COMM_WORLD,KaaS,ierr); CHKERRQ(ierr)
  call MatSetSizes(KaaS,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr); CHKERRQ(ierr)
  call MatSetType(KaaS,MATMPISBAIJ,ierr); CHKERRQ(ierr)
  call MatSetFromOptions(KaaS,ierr); CHKERRQ(ierr)
  call MatSetUp(KaaS,ierr); CHKERRQ(ierr)
  
  allocate(ranges(size+1))

  call MatGetOwnershipRange(KaaS, start, ende, ierr); CHKERRQ(ierr)
  
  ranges(rank+1) = start+1  ! start ist ein C-Feldindex und beginnt somit bei 0
  if (rank .eq. size-1) ranges(size+1) = ende+1
  
  do ii = 1,size
    call MPI_Bcast ( ranges(ii), 1, MPI_INTEGER, ii-1, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
  end do
  
  call MPI_Bcast ( ranges(size+1), 1, MPI_INTEGER, size-1, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  ENDE  - Bereite Matrizen vor'

  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  START - Bestimme Anzahl der Nullelemente'
  allocate(d_nzz(fesim%internals%dim_dof), o_nzz(fesim%internals%dim_dof))

  if (rank .eq. 0) then
  
    call SYSTEM_CLOCK(st,ctr,ctm)
    call estimate_number_nonzeros(fesim, virtNodes, ranges, size, d_nzz, o_nzz, .false.)
    call SYSTEM_CLOCK(en,ctr,ctm)
    write(*,*) '  Schaetzzeit: ', (en-st)/DBLE(ctr)
    
    
    do ii = 1,fesim%internals%dim_dof
      d_nzz(ii) = MIN(d_nzz(ii), fesim%internals%dim_dof-ii+1)
      o_nzz(ii) = MIN(o_nzz(ii), fesim%internals%dim_dof-ii+1)
    end do
  end if
  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  ENDE  - Bestimme Anzahl der Nullelemente'

  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  START - Lege Matrix an'
  call MPI_Bcast ( d_nzz(1:fesim%internals%dim_dof), fesim%internals%dim_dof, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
  call MPI_Bcast ( o_nzz(1:fesim%internals%dim_dof), fesim%internals%dim_dof, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
  
  bs = 1
  call MatMPISBAIJSetPreallocation(KaaS, bs, PETSC_NULL_INTEGER, d_nzz(ranges(rank+1):ranges(rank+2)-1), PETSC_NULL_INTEGER, o_nzz(ranges(rank+1):ranges(rank+2)-1), ierr); CHKERRQ(ierr)

  deallocate(d_nzz, o_nzz)
  deallocate(ranges)
  if (textoutput .eq. .true. .and. rank .eq. 0) write(*,*) '  ENDE  - Lege Matrix an'

  if (fesim%sol .eq. 2 .and. KgaaS .ne. PETSC_NULL_MAT) then
  
    call MatCreate(PETSC_COMM_WORLD,KgaaS,ierr); CHKERRQ(ierr)
    m = fesim%internals%dim_dof
    call MatSetSizes(KgaaS,PETSC_DECIDE,PETSC_DECIDE,m,m,ierr); CHKERRQ(ierr)
    call MatSetType(KgaaS,MATMPISBAIJ,ierr); CHKERRQ(ierr)
    call MatSetFromOptions(KgaaS,ierr); CHKERRQ(ierr)
    !call MatSetOption(KgaaS,MAT_IGNORE_LOWER_TRIANGULAR,PETSC_TRUE,ierr); CHKERRQ(ierr)
    call MatSetUp(KgaaS,ierr); CHKERRQ(ierr)
    
    allocate(ranges(size+1))

    call MatGetOwnershipRange(KgaaS, start, ende, ierr); CHKERRQ(ierr)
    
    ranges(rank+1) = start+1  ! start ist ein C-Feldindex und beginnt somit bei 0
    
    if (rank .eq. size-1) ranges(size+1) = ende+1
    
    do ii = 1,size
      call MPI_Bcast ( ranges(ii), 1, MPI_INTEGER, ii-1, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end do
    
    call MPI_Bcast ( ranges(size+1), 1, MPI_INTEGER, size-1, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    
    allocate(d_nzz(fesim%internals%dim_dof), o_nzz(fesim%internals%dim_dof))
    
    if (rank .eq. 0) then
      call estimate_number_nonzeros(fesim, virtNodes, ranges, size, d_nzz, o_nzz, .false.)
      do ii = 1,fesim%internals%dim_dof
        d_nzz(ii) = MIN(d_nzz(ii), fesim%internals%dim_dof)
        o_nzz(ii) = MIN(o_nzz(ii), fesim%internals%dim_dof)
      end do
    end if
    
    call MPI_Bcast ( d_nzz(1:fesim%internals%dim_dof), fesim%internals%dim_dof, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    call MPI_Bcast ( o_nzz(1:fesim%internals%dim_dof), fesim%internals%dim_dof, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    
    bs = 1
    call MatMPISBAIJSetPreallocation(KgaaS, bs, PETSC_NULL_INTEGER, d_nzz(ranges(rank+1):ranges(rank+2)-1), PETSC_NULL_INTEGER, o_nzz(ranges(rank+1):ranges(rank+2)-1), ierr); CHKERRQ(ierr)
    
    allocate(d_nzz_kgeo(fesim%internals%dim_dof), o_nzz_kgeo(fesim%internals%dim_dof))
    
    d_nzz_kgeo(1:fesim%internals%dim_dof) = d_nzz(1:fesim%internals%dim_dof)
    o_nzz_kgeo(1:fesim%internals%dim_dof) = o_nzz(1:fesim%internals%dim_dof)
    
    deallocate(d_nzz, o_nzz)
    deallocate(ranges)
    
  end if
  
  if (rank .eq. 0) then
    call clear_node_node_con(virtNodes,fesim%num_nodes)
    deallocate(virtNodes)
  end if
!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'init_matrices'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if
  
end subroutine init_matrices
