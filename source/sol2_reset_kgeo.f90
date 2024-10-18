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
!> $Id: sol2_reset_kgeo.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol2_reset_kgeo(dim_dof, KgaaS)
!
! Use
!
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
! Data types
!
  integer, intent(in) :: dim_dof
  integer        :: ii, size, rank
  
  integer, dimension(:), allocatable      :: ranges 

  
  PetscInt       :: m
  PetscInt       :: bs
  PetscErrorCode :: ierr
  PetscInt       :: start, ende

  Mat            :: KgaaS
  
  integer        :: err_code=0
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
  call MPI_Comm_size(PETSC_COMM_WORLD,size,ierr); CHKERRQ(ierr)
  call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
  
  call MatCreate(PETSC_COMM_WORLD,KgaaS,ierr); CHKERRQ(ierr)
  m = dim_dof
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
  
  bs = 1
  call MatMPISBAIJSetPreallocation(KgaaS, bs, PETSC_NULL_INTEGER, d_nzz_kgeo(ranges(rank+1):ranges(rank+2)-1), PETSC_NULL_INTEGER, o_nzz_kgeo(ranges(rank+1):ranges(rank+2)-1), ierr); CHKERRQ(ierr)

  

  deallocate(ranges)

!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol2_reset_kgeo'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if
  
end subroutine sol2_reset_kgeo
