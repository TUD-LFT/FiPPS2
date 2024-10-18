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
subroutine get_restored_eigenvektors(eigvals, eigenvectors_cond, eigenvectors, nconv, dim_dof_free, dim_dof, condens, rank, Kii, Kgii, Kib, Kgib)

#include "petsc/finclude/petscksp.h"
  use petscksp

  implicit none

  integer                                         :: dim_dof, dim_dof_free
  PetscMPIInt                                     :: rank
  PetscInt                                        :: nconv
  Mat                                             :: Kii, Kgii, Kib, Kgib
  
  double precision, dimension(nconv)              :: eigvals
  double precision, dimension(dim_dof_free,nconv) :: eigenvectors_cond
  double precision, dimension(dim_dof,nconv)      :: eigenvectors
  double precision, dimension(dim_dof-dim_dof_free) :: eigenvectors_free
  integer, dimension(dim_dof)                     :: condens
  Mat                                             :: Kib_lamda_KGib, Kii_lamda_KGii, invKii_lamda_KGii, CMat
  VecScatter                                      :: scatter
  Vec                                             :: qi, qb, qi_seq
  PetscScalar, pointer                            :: xx_v(:)
  integer                                         :: ii, mm, nn, i0, i1
  PetscErrorCode                                  :: ierr
  
  eigenvectors=0.d0
  do nn = 1, nconv
  
    write(*,*) 'Eigenvektor', nn
  
    write(*,*) 'Matrix Kgib kopieren'
    call MatConvert(Kgib, MATMPIAIJ, MAT_INITIAL_MATRIX, Kib_lamda_Kgib, ierr); CHKERRQ(ierr)
    write(*,*) 'Bilden von Kib-lamba*Kgib'
    call MatAYPX(Kib_lamda_Kgib, -eigvals(nn), Kib, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
    
    write(*,*) 'Matrix Kgii kopieren'
    call MatConvert(Kgii, MATMPIAIJ, MAT_INITIAL_MATRIX, Kii_lamda_Kgii, ierr); CHKERRQ(ierr)
    write(*,*) 'Bilden von Kii-lamba*Kgii'
    call MatAYPX(Kii_lamda_Kgii, -eigvals(nn), Kii, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
    
    write(*,*) 'Bilden von inv(Kii-lamba*Kgii)'
    call get_inverse3(Kii_lamda_Kgii, dim_dof-dim_dof_free, invKii_lamda_KGii, .FALSE.)
    
    write(*,*) 'Bilden von inv(Kii-lamba*Kgii)(Kib-lamba*Kgib)'
    call MatMatMult(invKii_lamda_KGii, Kib_lamda_Kgib, MAT_INITIAL_MATRIX, 1.d0, CMat, ierr); CHKERRQ(ierr)
    
    write(*,*) 'Bilden von -inv(Kii-lamba*Kgii)(Kib-lamba*Kgib)'
    call MatScale(CMat, -1.d0, ierr); CHKERRQ(ierr)
    
    write(*,*) 'Vektor qb erstellen'
    call VecCreate(PETSC_COMM_WORLD, qb, ierr); CHKERRQ(ierr)
    call VecSetSizes(qb,PETSC_DECIDE, dim_dof_free, ierr); CHKERRQ(ierr)
    call VecSetFromOptions(qb, ierr); CHKERRQ(ierr)
    
    do mm = 1, dim_dof_free
      call VecSetValue(qb, mm-1, eigenvectors_cond(mm, nn), INSERT_VALUES, ierr); CHKERRQ(ierr)
    end do
    
    call VecAssemblyBegin(qb,ierr); CHKERRQ(ierr);
    call VecAssemblyEnd(qb,ierr); CHKERRQ(ierr);
    
    call VecCreate(PETSC_COMM_WORLD, qi, ierr); CHKERRQ(ierr)
    call VecSetSizes(qi,PETSC_DECIDE, dim_dof-dim_dof_free, ierr); CHKERRQ(ierr)
    call VecSetFromOptions(qi, ierr); CHKERRQ(ierr)
    call VecAssemblyBegin(qi,ierr); CHKERRQ(ierr);
    call VecAssemblyEnd(qi,ierr); CHKERRQ(ierr);
    
    write(*,*) 'Bilden von -inv(Kii-lamba*Kgii)(Kib-lamba*Kgib)*qb'
    call MatMult(CMat, qb, qi, ierr); CHKERRQ(ierr)
    
    write(*,*) 'Vektor qi uebertragen'
    call VecScatterCreateToZero(qi, scatter, qi_seq, ierr); CHKERRQ(ierr)
   
    call VecScatterBegin(scatter, qi, qi_seq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    call VecScatterEnd  (scatter, qi, qi_seq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
    
    call VecScatterDestroy(scatter, ierr); CHKERRQ(ierr)
           
    if (rank .eq. 0) then

      call VecGetArrayF90(qi_seq, xx_v, ierr); CHKERRQ(ierr)
      do ii = 1, dim_dof-dim_dof_free
        eigenvectors_free(ii) = xx_v(ii)
      end do

      call VecRestoreArrayF90(qi_seq, xx_v, ierr); CHKERRQ(ierr)

    end if
           
    call VecDestroy(qi_seq, ierr); CHKERRQ(ierr)

    write(*,*) 'Einsortieren von qb und qi'
    i0 = 0
    i1 = 0
    do mm = 1, dim_dof
      if (condens(mm) .eq. 1) then
    	i1 = i1 + 1
    	eigenvectors(mm,nn) = eigenvectors_free(i1)
	!eigenvectors(mm,nn) = 0.d0
      else
    	i0 = i0 + 1
    	eigenvectors(mm,nn) = eigenvectors_cond(i0,nn)
      end if
    end do
    
    call MatDestroy(Kib_lamda_KGib, ierr); CHKERRQ(ierr)
    call MatDestroy(Kii_lamda_KGii, ierr); CHKERRQ(ierr)
    call MatDestroy(invKii_lamda_KGii, ierr); CHKERRQ(ierr)
    call MatDestroy(CMat, ierr); CHKERRQ(ierr)
    call VecDestroy(qb, ierr); CHKERRQ(ierr)
    call VecDestroy(qi, ierr); CHKERRQ(ierr)
    
  end do

end subroutine get_restored_eigenvektors
