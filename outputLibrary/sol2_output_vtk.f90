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
!> VTK-Ergebnisausgabe fuer lineare Stabilitaetsanalyse
!
!> @details
!
!> @author Florian Dexl, TU Dresden, WiMi, 30.09.2021
!
!> $Id: sol1_output_hymowi2.f90 434 2021-01-04 16:50:30Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 434 $
!> $Date: 2021-01-04 17:50:30 +0100 (Mo, 04. Jan 2021) $
!
! =================================================================================================
subroutine sol2_output_vtk(fesim, scloop, eigenvalues, eigenvectors, num_out_eigval, aeroConverged)
!
! use
!
  use fesimulation_typen
  use konstanten
  
  implicit none
  
  type(fe_simulation), intent(in)                         :: fesim
  double precision, dimension(num_out_eigval), intent(in) :: eigenvalues
  double precision, dimension(fesim%internals%dim_dof,num_out_eigval), intent(in) :: eigenvectors
  integer, intent(in)                                     :: num_out_eigval
  integer,intent(in)                                      :: scloop
  logical, intent(in)                                     :: aeroConverged
!
! Internal
!
  integer                                                 :: kk, ii, iii, pii, mm
  double precision, dimension(:),allocatable              :: U
  
  ! Routine

  allocate(U(fesim%num_dof))

  do kk = 1,num_out_eigval

    do ii=1, fesim%num_dof
       
      if (fesim%internals%num_dof_vec(ii) .GT. 0) then
        U(ii) = eigenvectors(fesim%internals%num_dof_vec(ii),kk)
      else if (fesim%internals%num_dof_vec(ii) .LT. 0) then
        U(ii) = 0.d0
        iii = fesim%internals%num_dof_vec(ii)
        do mm = 1,size(fesim%randbedingungen%mpcs(-iii)%independend,1)
          if (fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof .gt. 0) then !Freiheit ist durch SPC gesperrt, bereits aus dem Gleichungssystem ausgebaut und daher 0
            pii = fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(-iii)%independend(mm)%nid-1)*ndof+fesim%randbedingungen%mpcs(-iii)%independend(mm)%dof)
            U(ii) = U(ii) + eigenvectors(pii,kk)*fesim%randbedingungen%mpcs(-iii)%independend(mm)%fac
          end if
        end do
      else
        U(ii) = 0.d0
      end if
   
    end do
    
    do ii = 1,fesim%num_nodes
      if (fesim%knoten%nodes(ii)%cid .ne. 0) then
        U((ii-1)*ndof+1:(ii-1)*ndof+3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),U((ii-1)*ndof+1:(ii-1)*ndof+3))
        U((ii-1)*ndof+4:(ii-1)*ndof+6) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),U((ii-1)*ndof+4:(ii-1)*ndof+6))
      end if
    end do  
       
    call vtkouteig(fesim, kk, eigenvalues(kk), u)
    
  end do

  deallocate(U)

end subroutine sol2_output_vtk
