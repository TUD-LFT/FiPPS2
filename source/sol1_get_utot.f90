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
!> $Id: sol1_get_utot.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_get_Utot(fesim, Utot, UaS)

  use konstanten
  use fesimulation_typen
!
! =================================================================================================
!
! Include
!
#include "petsc/finclude/petscvec.h"
  use petscvec
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Data types
!
  type(fe_simulation)                   :: fesim
  integer                               :: ii, iii, pii, mm

  double precision, dimension(fesim%num_dof) :: Utot
  
  Vec                                   :: UaS
  
  PetscScalar, pointer                  :: xx_v(:)
  PetscErrorCode                        :: ierr
  
  integer                               :: err_code=0

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
  call VecGetArrayF90(UaS,xx_v,ierr); CHKERRQ(ierr)

  do ii=1, fesim%num_dof
    
    if (fesim%internals%num_dof_vec(ii) .GT. 0) then
      Utot(ii) = xx_v(fesim%internals%num_dof_vec(ii))
    else if (fesim%internals%num_dof_vec(ii) .LT. 0) then
      Utot(ii) = 0.d0
      iii = fesim%internals%num_dof_vec(ii)
      do mm = 1,size(fesim%randbedingungen%mpcs(-iii)%independend,1)
        if (fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof .gt. 0) then !Freiheit ist durch SPC gesperrt, bereits aus dem Gleichungssystem ausgebaut und daher 0
          pii = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-iii)%independend(mm)%nid-1)*ndof+fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof)
          Utot(ii) = Utot(ii) + xx_v(pii)*fesim%randbedingungen%mpcs(-iii)%independend(mm)%fac
        end if
      end do
    else
      Utot(ii) = 0.d0
    end if
    
  end do
  
  call VecRestoreArrayF90(UaS,xx_v,ierr); CHKERRQ(ierr)

  do ii = 1,fesim%num_nodes
    if (fesim%knoten%nodes(ii)%cid .ne. 0) then
      Utot((ii-1)*ndof+1:(ii-1)*ndof+3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat,Utot((ii-1)*ndof+1:(ii-1)*ndof+3))
      Utot((ii-1)*ndof+4:(ii-1)*ndof+6) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat,Utot((ii-1)*ndof+4:(ii-1)*ndof+6))
    end if
  end do

!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'sol1_get_Utot'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if


end subroutine sol1_get_Utot
