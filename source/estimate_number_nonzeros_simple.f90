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
!> $Id: estimate_number_nonzeros_simple.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine estimate_number_nonzeros_simple( fesim, virtNodes, ranges, mpi_size, d_nzz, o_nzz, geometry)
! =================================================================================================
!
! use
!
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
!
! =================================================================================================
!
! Data types
!
! Input
!
  type(fe_simulation)                                        :: fesim
  type(virtualNode), dimension(fesim%num_nodes),target, intent(in) :: virtNodes !< node pointer list
  integer, intent(in)                                        :: mpi_size        !< number of MPI ranks
  integer, intent(in), dimension(mpi_size+1)                 :: ranges          !< array which contains the index range of the global stiffness matrix for each mpi rank
  logical, intent(in)                                        :: geometry
!
! Output
!
  PetscInt, intent(out), dimension(fesim%internals%dim_dof)  :: d_nzz           !< array containing the number of block nonzeros in the various block rows in the upper triangular and diagonal part of the in diagonal portion of the local --> d_nzz
  PetscInt, intent(out), dimension(fesim%internals%dim_dof)  :: o_nzz           !< array containing the number of nonzeros in the various block rows of the off-diagonal portion of the local submatrix --> o_nzz
!
! Internal
!
  integer                                                    :: ii, iii, jjj, node_dof, range
  integer                                                    :: actNodeID, mm, nn
  type(virtualNode), Pointer                                 :: actNode  
  
  integer                                                    :: err_code=0
  integer, allocatable, dimension(:)                         :: iiVek, jjVek
!
! =================================================================================================
!
! Initialisation
!
  d_nzz = 0
  o_nzz = 0
!
! =================================================================================================
!
! Calculation
!

  do actNodeID = 1, fesim%num_nodes
    do node_dof = 1,ndof
      iii = fesim%internals%num_dof_vec((actNodeID-1)*ndof+node_dof)
      if (iii .eq. 0) cycle
      
      if (iii .GT. 0) then
        allocate(iiVek(1))
        iiVek(1) = iii
      else
        allocate(iiVek(size(fesim%randbedingungen%mpcs(-iii)%independend,1)))
        do mm = 1,size(iiVek,1)
          if (fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof .gt. 0) then
            iiVek(mm) = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-iii)%independend(mm)%nid-1)*ndof+fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof)
          else
            iiVek(mm) = -10
          end if
        end do
      end if
      
      DO mm = 1,size(iiVek,1)
        iii = iiVek(mm)
        if (iii .eq. -10) cycle
        range = 1
        do while(iii .ge. ranges(range+1))
          range = range + 1
        end do
        actNode => virtNodes(actNodeID)
        DO WHILE(ASSOCIATED(actNode))
          do ii = 1,ndof
            if (geometry .EQ. .TRUE.) then
              if (node_dof .GT. actNode%numDofGeo .OR. ii .GT. actNode%numDofGeo) cycle
            end if
            jjj = fesim%internals%num_dof_vec((actNode%nodeID-1)*ndof+ii)
            if (jjj  .eq. 0) cycle
            if (jjj .GT. 0) then
              allocate(jjVek(1))
              jjVek(1) = jjj
            else
              allocate(jjVek(size(fesim%randbedingungen%mpcs(-jjj)%independend,1)))
              do nn = 1,size(jjVek,1)
                if (fesim%randbedingungen%mpcs(-jjj)%independend(nn)%dof .gt. 0) then
                  jjVek(nn) = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-jjj)%independend(nn)%nid-1)*ndof+fesim%randbedingungen%mpcs(-jjj)%independend(nn)%dof)
                else
                  jjVek(nn) = -10
                end if
              end do
            end if
   
            DO nn = 1, size(jjVek,1)
              jjj = jjVek(nn)
              if (jjj .eq. -10) cycle
              if (jjj .lt. iii) cycle
              if (jjj .ge. ranges(range+1)) then
                o_nzz(iii) = o_nzz(iii) + 1
              else
                d_nzz(iii) = d_nzz(iii) + 1
              end if
            END DO
            deallocate(jjVek)
          end do
          actNode => actNode%nextNode
        END DO
      END DO
      deallocate(iiVek)
    end do
  end do
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'estimate_number_nonzeros_simple'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine estimate_number_nonzeros_simple
