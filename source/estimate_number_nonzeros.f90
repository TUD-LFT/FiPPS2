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
!> Subroutine estimates the number of nonzero-enries per row in global 
!> stiffness matrix for worst-case (element matrices fully occupied)
!> This number is needed for sufficient memory allocation in PETSc.
!> To get a better estimation the number of nonzero elements is determined
!> for portions of each row. The number of elements inside one portion is
!> blocksize.
!
!> @details
!> Subroutine estimates the number of nonzero-enries per row in global 
!> stiffness matrix for worst-case (element matrices fully occupied)
!> This number is needed for sufficient memory allocation in PETSc.
!> For that purpose the virtNodes pointer list.
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter 22.06.2010
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 20.09.2010
!
!> $Id: estimate_number_nonzeros.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine estimate_number_nonzeros( fesim, virtNodes, ranges, mpi_size, d_nzz, o_nzz, geometry)
! =================================================================================================
!
! use
!
  use konstanten
  use fesimulation_typen
  use pre_assemble_types
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscmat.h"
  use petscmat
!
! =================================================================================================
!
  implicit none!
!
!
  INTERFACE 
    SUBROUTINE ESTIMATE_NUMBER_NONZEROS_SIMPLE(FESIM,VIRTNODES,RANGES, &
    & MPI_SIZE,D_NZZ,O_NZZ,GEOMETRY)
      USE FESIMULATION_TYPEN
      USE PRE_ASSEMBLE_TYPES
      TYPE(FE_SIMULATION) :: FESIM                                      
      INTEGER(KIND=4), INTENT(IN) :: MPI_SIZE                         !< number of MPI ranks
      TYPE (VIRTUALNODE) ,TARGET, INTENT(IN) :: VIRTNODES(      &
              &FESIM%NUM_NODES)                                       !< node pointer list 
      INTEGER(KIND=4), INTENT(IN) :: RANGES(MPI_SIZE+1)               !< array containing the number of block nonzeros in the various block rows in the upper triangular and diagonal part of the in diagonal portion of the local
      INTEGER(KIND=4), INTENT(OUT) :: D_NZZ(FESIM%INTERNALS%DIM_DOF)  !< array containing the number of block nonzeros in the various block rows in the upper triangular and diagonal part of the in diagonal portion of the local
      INTEGER(KIND=4), INTENT(OUT) :: O_NZZ(FESIM%INTERNALS%DIM_DOF)  !< array containing the number of nonzeros in the various block rows of the off-diagonal portion of the local submatrix 
      LOGICAL(KIND=4), INTENT(IN) :: GEOMETRY
    END SUBROUTINE ESTIMATE_NUMBER_NONZEROS_SIMPLE
  END INTERFACE  
  INTERFACE 
    SUBROUTINE ESTIMATE_NUMBER_NONZEROS_BLOCK(FESIM,VIRTNODES,RANGES, &
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
    END SUBROUTINE ESTIMATE_NUMBER_NONZEROS_BLOCK
  END INTERFACE 
!=================================================================================================
!
! Data types
!
! Input
!
  type(fe_simulation)                                        :: fesim
  type(virtualNode), dimension(fesim%num_nodes),target, intent(in) :: virtNodes
  integer, intent(in)                                        :: mpi_size
  integer, intent(in), dimension(mpi_size+1)                 :: ranges
  logical, intent(in)                                        :: geometry
!
! Output
!
  PetscInt, intent(out), dimension(fesim%internals%dim_dof)  :: d_nzz, o_nzz
!
! Internal
!
  integer                                                    :: err_code=0
!
! =================================================================================================
!
! Calculation
!
  if (fesim%blocksize .GE. fesim%internals%dim_dof) then
    if (textoutput .eq. .true.) write(*,*) '   einfache Abschaetzung'
    call estimate_number_nonzeros_simple( fesim, virtNodes, ranges, mpi_size, d_nzz, o_nzz, geometry)
  else
    if (textoutput .eq. .true.) write(*,*) '   Blockabschaetzung'
    if (textoutput .eq. .true.) write(*,*) '     blocksize: ', fesim%blocksize
    call estimate_number_nonzeros_block( fesim, virtNodes, ranges, mpi_size, d_nzz, o_nzz, geometry)
  end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'estimate_number_nonzeros'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine estimate_number_nonzeros
