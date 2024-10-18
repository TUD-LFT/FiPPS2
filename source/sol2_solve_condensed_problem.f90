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
!> subroutine for computing load-vectors entrys resulting from single forces
!
!> @details
!> subroutine calculates the force contribution to the loadvector
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit, 08.12.2009
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit, 09.10.2015
!>                                  überflüssige Schleife über nodes(:) entfernt,
!>                                  zur Einsparung von Rechenzeit
!
!> $Id: load_force.f90 411 2019-10-16 13:38:42Z DOM+ahauffe $
!> $Author: DOM+ahauffe $
!> $Revision: 411 $
!> $Date: 2019-10-16 15:38:42 +0200 (Mi, 16. Okt 2019) $
!
! =================================================================================================
subroutine sol2_solve_condensed_problem (rank,fesim,ksp,KaaS,KgaaS,scloop)
!
! Use
!
#include "slepc/finclude/slepceps.h"
  use slepceps
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
  type(fe_simulation)                           :: fesim
  
  PetscMPIInt                                   :: rank
  
  KSP                                           :: ksp

  Mat                                           :: KaaS, KaaSMPIAIJ
  Mat                                           :: KgaaS, KgaaSMPIAIJ

  integer                                       :: err_code=0
  integer, parameter                            :: fileunit = 83
  
  !double precision, dimension(:,:), allocatable :: condensed_freedoms

  IS                                            :: is0, is1
  
  integer, dimension(:), allocatable            :: condens 
  integer                                       :: i0, i1, mm, nn
  PetscInt, dimension(:), allocatable           :: index0, index1
  integer                                       :: numCondNodes, nodeID, numCondDOFs
  Mat                                           :: KSchur, KgSchur
  
  Mat                                           :: SubKbb, SubKib, SubKbi, SubKii, SubInvKiiSparse
  Mat                                           :: SubKgbb, SubKgib, SubKgbi, SubKgii, SubInvKgiiSparse
  Mat                                           :: T, SubKbiT, SubKgbiT
  PetscReal                                     :: fill
  
  EPS                                           :: eps
  PetscInt                                      :: nconv
  double precision, dimension(:), allocatable   :: eigvals
  double precision, dimension(:,:), allocatable :: eigenvectors, eigenvectors_cond
  PetscErrorCode                                :: ierr
  integer                                       :: scloop
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation
     
  open(unit=fileunit, file='nodestocondens.fipps', status='old', action='read')
      
  read(fileunit,'(I10)') numCondNodes

  allocate(condens(fesim%internals%dim_dof))

  condens = 0

  numCondDOFs = 0

!  allocate(condensed_freedoms(fesim%num_nodes,6))
!  condensed_freedoms = 0.d0

  do mm = 1, numCondNodes
    read(fileunit,'(I10)') nodeID
    do nn = 1,6
      if (fesim%internals%num_dof_vec((nodeID-1)*ndof+nn) .NE. 0 .AND. condens(fesim%internals%num_dof_vec((nodeID-1)*ndof+nn)) .NE. 1) then
        condens(fesim%internals%num_dof_vec((nodeID-1)*ndof+nn)) = 1
!        condensed_freedoms(nodeID,nn) = 1.0
        numCondDOFs = numCondDOFs + 1
      end if
    end do
  end do

!  do nn = 1,6
!    write(tempChar,'(I1)') nn
!    call vtkoutpointscalar(condensed_freedoms(:,nn), num_nodes, 'dof' // tempChar // '                          ')
!  end do

!  deallocate(condensed_freedoms)

  close(fileunit)
  
  write(*,*) 'dim_dof-numCondDOFs 1', fesim%internals%dim_dof-numCondDOFs
  allocate(index0(numCondDOFs), index1(fesim%internals%dim_dof-numCondDOFs))

  i0 = 0; i1 = 0

  do mm = 1,fesim%internals%dim_dof 
    if (condens(mm) .EQ. 1) then
      i0 = i0 + 1
      index0(i0) = mm-1
    else
      i1 = i1 + 1
      index1(i1) = mm-1
    end if
  end do
  
  call ISCreateGeneral(PETSC_COMM_WORLD,numCondDOFs,index0,PETSC_COPY_VALUES,is0,ierr); CHKERRQ(ierr) ! zu kondensierende DOF (Interior)
  call ISCreateGeneral(PETSC_COMM_WORLD,fesim%internals%dim_dof-numCondDOFs,index1,PETSC_COPY_VALUES,is1,ierr); CHKERRQ(ierr) ! freie DOF
  
  deallocate(index0, index1)

  write(*,*) 'Steifigkeitsmatrix kondensieren'
  ! Erzeugen der Submatrizen
  ! | Kbb Kbi |   | ub |   | fb |
  ! |	      | * |    | = |	|
  ! | Kib Kii |   | ui |   | fi |
  WRITE(*,*) 'Erzeugen der Submatrizen'
  
  call MatConvert(KaaS, MATMPIAIJ, MAT_INITIAL_MATRIX, KaaSMPIAIJ, ierr); CHKERRQ(ierr)
  
  call MatCreateSubMatrix(KaaSMPIAIJ, is1, is1, MAT_INITIAL_MATRIX, SubKbb, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(KaaSMPIAIJ, is0, is0, MAT_INITIAL_MATRIX, SubKii, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(KaaSMPIAIJ, is1, is0, MAT_INITIAL_MATRIX, SubKbi, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(KaaSMPIAIJ, is0, is1, MAT_INITIAL_MATRIX, SubKib, ierr); CHKERRQ(ierr)
  
  call MatDestroy(KaaSMPIAIJ, ierr); CHKERRQ(ierr)
  
  ! Berechnen von inv(Kii) als Densematrix
  WRITE(*,*) 'Berechnen von inv(Kii) als Densematrix'
  call get_inverse3(SubKii, numCondDOFs, SubInvKiiSparse, .false.)
  
  fill = 1.d0
  
  ! Berechnen von Kbi * Inv (Kii)
  WRITE(*,*) 'Berechnen von Kbi * Inv (Kii)'
  
  call MatMatMult(SubKbi, SubInvKiiSparse, MAT_INITIAL_MATRIX, fill, T, ierr); CHKERRQ(ierr)
  call MatDestroy(SubInvKiiSparse, ierr); CHKERRQ(ierr)

  ! Berechnen von Kbi * Inv (Kii) * Kbi^T
  WRITE(*,*) 'Berechnen von Kbi * Inv (Kii) * Kbi^T'

  call MatTranspose(SubKbi, MAT_INITIAL_MATRIX, SubKbiT, ierr); CHKERRQ(ierr)
  call MatMatMult(T, SubKbiT, MAT_INITIAL_MATRIX, fill, KSchur, ierr); CHKERRQ(ierr)
  call MatDestroy(T, ierr); CHKERRQ(ierr)
  call MatDestroy(SubKbiT, ierr); CHKERRQ(ierr)

  ! Berechnen von  -Kbb + Kbi * Inv (Kii) * Kbi^T
  WRITE(*,*) 'Berechnen von  -Kbb + Kbi * Inv (Kii) * Kbi^T'

  call MatAXPY(KSchur, -1.d0, SubKbb, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)

  ! Berechnen von Kbb - Kbi * Inv (Kii) * Kbi^T 
  WRITE(*,*) 'Berechnen von Kbb - Kbi * Inv (Kii) * Kbi^T'

  call MatScale(KSchur, -1.d0, ierr); CHKERRQ(ierr)  

  write(*,*) 'geometrische Steifigkeitsmatrix kondensieren'
  ! Erzeugen der Submatrizen
  ! | Kgbb Kgbi |   | ub |   | fb |
  ! |	      | * |    | = |	|
  ! | Kgib Kgii |   | ui |   | fi |
  WRITE(*,*) 'Erzeugen der Submatrizen'
  
  call MatConvert(KgaaS, MATMPIAIJ, MAT_INITIAL_MATRIX, KgaaSMPIAIJ, ierr); CHKERRQ(ierr)
  
  call MatCreateSubMatrix(KgaaSMPIAIJ, is1, is1, MAT_INITIAL_MATRIX, SubKgbb, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(KgaaSMPIAIJ, is0, is0, MAT_INITIAL_MATRIX, SubKgii, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(KgaaSMPIAIJ, is1, is0, MAT_INITIAL_MATRIX, SubKgbi, ierr); CHKERRQ(ierr)
  call MatCreateSubMatrix(KgaaSMPIAIJ, is0, is1, MAT_INITIAL_MATRIX, SubKgib, ierr); CHKERRQ(ierr)
  
  call MatDestroy(KgaaSMPIAIJ, ierr); CHKERRQ(ierr)

  call get_inverse3(SubKgii, numCondDOFs, SubInvKgiiSparse, .TRUE.)
  
  fill = 1.d0
  
  ! Berechnen von Kgbi * Inv (Kgii)
  WRITE(*,*) 'Berechnen von Kgbi * Inv (Kgii)'
  
  call MatMatMult(SubKgbi, SubInvKgiiSparse, MAT_INITIAL_MATRIX, fill, T, ierr); CHKERRQ(ierr)
  call MatDestroy(SubInvKgiiSparse, ierr); CHKERRQ(ierr)

  ! Berechnen von Kgbi * Inv (Kgii) * Kgbi^T
  WRITE(*,*) 'Berechnen von Kgbi * Inv (Kgii) * Kgbi^T'

  call MatTranspose(SubKgbi, MAT_INITIAL_MATRIX, SubKgbiT, ierr); CHKERRQ(ierr)
  call MatMatMult(T, SubKgbiT, MAT_INITIAL_MATRIX, fill, KgSchur, ierr); CHKERRQ(ierr)
  call MatDestroy(T, ierr); CHKERRQ(ierr)
  call MatDestroy(SubKgbiT, ierr); CHKERRQ(ierr)

  ! Berechnen von  -Kgbb + Kgbi * Inv (Kgii) * Kgbi^T
  WRITE(*,*) 'Berechnen von  -Kgbb + Kgbi * Inv (Kgii) * Kgbi^T'

  call MatAXPY(KgSchur, -1.d0, SubKgbb, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)

  ! Berechnen von Kgbb - Kgbi * Inv (Kgii) * Kgbi^T 
  WRITE(*,*) 'Berechnen von Kgbb - Kgbi * Inv (Kgii) * Kgbi^T'

  call MatScale(KgSchur, -1.d0, ierr); CHKERRQ(ierr)
  
  
  
  
  
  write(*,*) 'Kg fertig, LC', scloop
  call sol2_solve(eps, ksp, KSchur, KgSchur, fesim%numEigVal)
  write(*,*) 'Eigenwertproblem geloest, LC', scloop 
  
  call MatDestroy(KSchur,ierr); CHKERRQ(ierr)
  call MatDestroy(KgSchur,ierr); CHKERRQ(ierr)  

  call EPSGetConverged(eps,nconv,ierr); CHKERRQ(ierr) 
  
  write(*,*) 'Anzahl der konvergierten Eigenvektoren: ', nconv
  
  nconv = fesim%numEigVal

  if (rank == 0) then  
    allocate(eigvals(nconv))
    allocate(eigenvectors_cond(fesim%internals%dim_dof-numCondDOFs,nconv))
  else
    allocate(eigvals(1))
    allocate(eigenvectors_cond(1,1))
  end if
  
  write(*,*) 'dim_dof-numCondDOFs 2', fesim%internals%dim_dof-numCondDOFs
  
  call sol2_extract_eigenvalues(eps, eigvals, eigenvectors_cond, nconv, fesim%internals%dim_dof-numCondDOFs, rank)
  call EPSDestroy(eps,ierr); CHKERRQ(ierr)
  
  allocate(eigenvectors(fesim%internals%dim_dof,nconv))
  
  call get_restored_eigenvektors(eigvals, eigenvectors_cond, eigenvectors, nconv, fesim%internals%dim_dof-numCondDOFs, fesim%internals%dim_dof, condens, rank, SubKii, SubKgii, SubKib, SubKgib)
  deallocate(condens)

  !if (rank == 0) call sol2_output_vtk(fesim, 1, eigvals, eigenvectors, nconv, .False.)
  
  if (rank == 0) call sol2_output_control(fesim, scloop, eigvals, eigenvectors, nconv, .False.)

  deallocate(eigvals)
  deallocate(eigenvectors_cond)
  deallocate(eigenvectors)
  
  call MatDestroy(SubKbb, ierr); CHKERRQ(ierr)  
  call MatDestroy(SubKii, ierr); CHKERRQ(ierr)  
  call MatDestroy(SubKbi, ierr); CHKERRQ(ierr)  
  call MatDestroy(SubKib, ierr); CHKERRQ(ierr)
  call MatDestroy(SubKgbb, ierr); CHKERRQ(ierr)  
  call MatDestroy(SubKgii, ierr); CHKERRQ(ierr)  
  call MatDestroy(SubKgbi, ierr); CHKERRQ(ierr)  
  call MatDestroy(SubKgib, ierr); CHKERRQ(ierr)

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol2_prepare_condensed_problem'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine sol2_solve_condensed_problem
