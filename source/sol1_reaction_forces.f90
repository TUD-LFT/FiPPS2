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
subroutine sol1_reaction_forces(fesim, Utot, Ftot, rank)

#include "petsc/finclude/petscmat.h"
    use petscmat
    use konstanten
    use fesimulation_typen

    implicit none
    
    type(fe_simulation)                                     :: fesim
    double precision, dimension(:), allocatable, intent(in) :: Utot
    double precision, dimension(:), allocatable, intent(in) :: Ftot
    PetscMPIInt, intent(in)                                 :: rank
    
    integer, dimension(:), allocatable                      :: temp_num_dof_vec
    integer                                                 :: temp_dim_dof
    double precision, dimension(:), allocatable             :: Utotlocal
    
    double precision, dimension(fesim%num_dof)              :: freact
    
    Mat                                                     :: KaaS         ! Steifigkeitsmatrix ohne gesperrte Freiheiten als Sparse-Matrix
    PetscErrorCode                                          :: ierr
    Vec                                                     :: Uvec
    Vec                                                     :: Fvec,Fvecseq
    PetscScalar                                             :: v
    PetscInt                                                :: pii
    VecScatter                                              :: scatter
    PetscScalar, pointer                                    :: xx_v(:)
    
    integer                                                 :: ii
    integer                                                 :: err_code=0
    
    !
    ! Sichern und Umspeichern des num_dof_vecs
    !
    ! muss neu berechnet werden, da keine Randbedinungen mehr vorhanden sind
    ! muss aber zurückgeschrieben werden, um den weiteren Programmverlauf nicht zu 
    ! zerstören, wie weitere Lastfälle oder Beulrechnung
    !
    
    if (rank == 0) then
    
      allocate ( temp_num_dof_vec(fesim%num_dof) )
      
      temp_num_dof_vec = fesim%internals%num_dof_vec
      
      do ii=1,fesim%num_dof
        fesim%internals%num_dof_vec(ii) = ii
      end do
      
      temp_dim_dof = fesim%internals%dim_dof
      
      fesim%internals%dim_dof = fesim%num_dof
    else
      deallocate ( fesim%internals%num_dof_vec )
    end if
    
    call bcast_internals(fesim%internals,fesim%num_dof)
    
    ! PETSC Verschiebungsvektor anlegen
    
    ! Anlegen von Uvec
    call VecCreate(PETSC_COMM_WORLD,Uvec,ierr); CHKERRQ(ierr)
    call VecSetSizes(Uvec,PETSC_DECIDE,fesim%num_dof,ierr); CHKERRQ(ierr)
    call VecSetFromOptions(Uvec,ierr); CHKERRQ(ierr)
    
    if (rank == 0) then
      
      allocate( Utotlocal(fesim%num_dof) )
    
      do ii = 1,fesim%num_nodes
        if (fesim%knoten%nodes(ii)%cid .eq. 0) then
          Utotlocal((ii-1)*ndof+1:(ii-1)*ndof+6) = Utot((ii-1)*ndof+1:(ii-1)*ndof+6)
        else
          Utotlocal((ii-1)*ndof+1:(ii-1)*ndof+3) = matmul(transpose(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat),Utot((ii-1)*ndof+1:(ii-1)*ndof+3))
          Utotlocal((ii-1)*ndof+4:(ii-1)*ndof+6) = matmul(transpose(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat),Utot((ii-1)*ndof+4:(ii-1)*ndof+6))
        end if
      end do
    
      do ii = 1, fesim%num_dof
        if (Utotlocal(ii) .NE. 0.d0) then
          v = Utotlocal(ii)
          pii = ii-1
          call VecSetValue(Uvec, pii, v, INSERT_VALUES, ierr); CHKERRQ(ierr)
        end if
      end do
      
      deallocate( Utotlocal ) 
      
    end if
    
    call VecAssemblyBegin(Uvec,ierr); CHKERRA(ierr);
    call VecAssemblyEnd(Uvec,ierr); CHKERRA(ierr);
    
    call VecDuplicate(Uvec,Fvec,ierr); CHKERRA(ierr);
    
    ! Steifigkeitsmatrix erstellen und aufbauen
    
    call init_matrices(fesim, KaaS, PETSC_NULL_MAT, rank)
    
    call sol1_get_stiffness_matrices(fesim, KaaS)
    
    call MatAssemblyBegin(KaaS,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
    call MatAssemblyEnd  (KaaS,MAT_FINAL_ASSEMBLY,ierr); CHKERRA(ierr)
    
    ! Lösungsvektor berechnen
    
    call MatMult(KaaS,Uvec,Fvec,ierr); CHKERRA(ierr)
    
    call MatDestroy(KaaS,ierr); CHKERRA(ierr)
    call VecDestroy(Uvec,ierr); CHKERRA(ierr)
    
    ! Lösungsvektor umschreiben
    
    call VecDuplicate(Fvec,Fvecseq,ierr); CHKERRA(ierr);
    
    call VecScatterCreateToZero(Fvec, scatter, Fvecseq, ierr); CHKERRA(ierr)
    
    call VecScatterBegin(scatter, Fvec, Fvecseq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRA(ierr)
    call VecScatterEnd  (scatter, Fvec, Fvecseq, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRA(ierr)
    
    call VecScatterDestroy(scatter, ierr); CHKERRA(ierr)
    
    if (rank == 0) then
    
      call VecGetArrayF90(Fvecseq,xx_v,ierr); CHKERRQ(ierr)
    
      if (fesim%globalReactForce .eq. .true.) then
        do ii = 1,fesim%num_nodes
          if (fesim%knoten%nodes(ii)%cid .eq. 0) then
            freact((ii-1)*ndof+1:(ii-1)*ndof+6) = xx_v((ii-1)*ndof+1:(ii-1)*ndof+6) !- ftot((ii-1)*ndof+1:(ii-1)*ndof+6)
          else
            ! wenn ein Knotenkoordinatensystem gegeben ist, dann die Verschiebungen aus dem lokalen ins globale Koordinatensystem drehen
            freact((ii-1)*ndof+1:(ii-1)*ndof+3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat,xx_v((ii-1)*ndof+1:(ii-1)*ndof+3)) !- ftot((ii-1)*ndof+1:(ii-1)*ndof+3)
            freact((ii-1)*ndof+4:(ii-1)*ndof+6) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat,xx_v((ii-1)*ndof+4:(ii-1)*ndof+6)) !- ftot((ii-1)*ndof+4:(ii-1)*ndof+6)
          end if
        end do
      else
        do ii = 1,fesim%num_nodes
          if (fesim%knoten%nodes(ii)%cid .eq. 0) then
            freact((ii-1)*ndof+1:(ii-1)*ndof+6) = xx_v((ii-1)*ndof+1:(ii-1)*ndof+6) - ftot((ii-1)*ndof+1:(ii-1)*ndof+6)
          else
            ! wenn ein Knotenkoordinatensystem gegeben ist, dann die aufgebrachten Kräfte aus dem globalen ins lokale Koordinatensystem drehen
            freact((ii-1)*ndof+1:(ii-1)*ndof+3) = xx_v((ii-1)*ndof+1:(ii-1)*ndof+3) !- matmul(transpose(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat),ftot((ii-1)*ndof+1:(ii-1)*ndof+3))
            freact((ii-1)*ndof+4:(ii-1)*ndof+6) = xx_v((ii-1)*ndof+4:(ii-1)*ndof+6) !- matmul(transpose(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat),ftot((ii-1)*ndof+4:(ii-1)*ndof+6))
          end if
        end do
      end if
      
    ! VTK-Ausgabe
      
      if (fesim%globalReactForce .eq. .true.) then
        call vtkoutstatic(fesim, freact, 'Reaktionskraefte_global', .false.)
        call vtkoutstatic(fesim, freact, 'Reaktionsmomente_global', .true.)
      else
        call vtkoutstatic(fesim, freact, 'Reaktionskraefte_lokal', .false.)
        call vtkoutstatic(fesim, freact, 'Reaktionsmomente_lokal', .true.)
      end if
    
    end if
    
    ! Aufräumen
    
    call VecDestroy(Fvec,ierr); CHKERRA(ierr)
    call VecDestroy(Fvecseq,ierr); CHKERRA(ierr)
    
    if (rank == 0) then
    
      fesim%internals%num_dof_vec = temp_num_dof_vec
      
      deallocate( temp_num_dof_vec )
      
      fesim%internals%dim_dof = temp_dim_dof
      
    else
      deallocate ( fesim%internals%num_dof_vec )
    end if
    
    call bcast_internals(fesim%internals,fesim%num_dof)
 
!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol1_reaction_forces'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine sol1_reaction_forces
