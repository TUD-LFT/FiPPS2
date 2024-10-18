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
!> control subroutine for inserting homogeneous boundary conditions
!
!> @details
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter 08.12.2009
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 11.10.2012
!
!> $Id: boundary_varload_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine boundary_varload_control (fesim,scloop)
! =================================================================================================
!
!       Header:         control subroutine for inserting homogeneous boundary conditions
!
!       Content:        
!
!       Input:          
!
!       Output:         
!
!       Internal:       
!
!       Calls:          
!
!       Called by:      
!
!       Author:         Martin Rädel                    08.12.2009
!                       TU Dresden, Diplomarbeit
!                       Andreas Hauffe                  11.10.2012
!                       TU Dresden, WM
!
!       Revision:       
!
! =================================================================================================
!
use konstanten
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!

!
! Output
!
!
! Input + Output
!
type(fe_simulation)                     :: fesim
integer, intent(in)                     :: scloop
!
! inner
!
integer                                 :: oldCpSet
integer                                 :: actDOF

integer                                 :: gg,ii
integer                                 :: err_code=0

integer                                 :: changeInd
logical, dimension(:), allocatable      :: fixedU
!
! =================================================================================================
!
! Calculation
!

allocate(fixedU(fesim%num_dof))

! number degrees of freedom sequencially

do ii=1,fesim%num_dof
   fesim%internals%num_dof_vec(ii) = ii
end do

fixedU = .FALSE.

! Set DOF for homogeneous bc to fixed
call boundary_spc_control (fesim,fixedU,scloop)

! Set DOF for couplings to fixed
call boundary_coupling_control (fesim,fixedU)

if (textoutput .eq. .true.) write(*,*) '  START - Integriere Platzhalter fuer MPCs'
! Set dependend DOF for MPC to fixed
call boundary_mpc_control (fesim,fixedU,scloop)
if (textoutput .eq. .true.) write(*,*) '  ENDE  - Integriere Platzhalter fuer MPCs'

if (textoutput .eq. .true.) write(*,*) '  START - Sortiere Freiheiten'
! Sort for free and fixed dofs
changeInd = fesim%num_dof
do ii = 1, fesim%num_dof

   if (changeInd == ii) exit
  
   if (fixedU(ii) == .TRUE.) then
  
      do while (fixedU(changeInd) == .TRUE.)
         changeInd = changeInd - 1
      end do
      
      if (changeInd .LT. ii) exit
      
      fixedU(ii)        = .FALSE.
      fixedU(changeInd) = .TRUE.
      fesim%internals%num_dof_vec(ii) = changeInd
      fesim%internals%num_dof_vec(changeInd) = ii
      changeInd = changeInd - 1
      
   end if
end do
if (textoutput .eq. .true.) write(*,*) '  ENDE  - Sortiere Freiheiten'


! get dimension of free dofs in equation system

fesim%internals%dim_dof = 0
do ii = 1, fesim%num_dof
  if (fixedU(ii) == .TRUE.) then
    exit
  else
    fesim%internals%dim_dof = fesim%internals%dim_dof + 1
  end if
end do

do ii = 1, fesim%num_dof
  if (fesim%internals%num_dof_vec(ii) .GT. fesim%internals%dim_dof) fesim%internals%num_dof_vec(ii) = 0
end do

if (textoutput .eq. .true.) write(*,*) '  START - Integriere Couplings'
! Add Coupling independend DOF-IDs to the dependend DOF for assembling
if (fesim%is_coupling .EQ. .true.) then
  oldCPSet = fesim%randbedingungen%couplings(1)%cpsid
  actDOF   = fesim%internals%num_dof_vec((fesim%randbedingungen%couplings(1)%nid-1)*ndof+fesim%randbedingungen%couplings(1)%dof)
  do ii = 2, size(fesim%randbedingungen%couplings,1)
    if (fesim%randbedingungen%couplings(ii)%cpsid .eq. oldCPSet) then
      fesim%internals%num_dof_vec((fesim%randbedingungen%couplings(ii)%nid-1)*ndof+fesim%randbedingungen%couplings(ii)%dof) = actDOF
    else
      oldCPSet = fesim%randbedingungen%couplings(ii)%cpsid
      actDOF   = fesim%internals%num_dof_vec((fesim%randbedingungen%couplings(ii)%nid-1)*ndof+fesim%randbedingungen%couplings(ii)%dof)
    end if
  end do
end if
if (textoutput .eq. .true.) write(*,*) '  ENDE  - Integriere Couplings'

if (textoutput .eq. .true.) write(*,*) '  START - Integriere MPCs'
! Add MPC-IDs to the dependend DOFs for assembling
! Loop over MPCs
if (fesim%is_mpcadd == .true.) then
    do gg=1,size(fesim%randbedingungen%mpcadds,1)
        if (fesim%lasten%subcases(scloop)%mpcaddid == fesim%randbedingungen%mpcadds(gg)%scid .or. &
          & fesim%randbedingungen%mpcadds(gg)%scid == -10) then                                             ! interne MPC durch Kontakte
            
            ! Multi point constraints
            
            if (fesim%is_mpc == .true.) then
                !do ii=1,size(mpcs,1)
                ii = fesim%randbedingungen%mpcadds(gg)%sid
                if (fesim%randbedingungen%mpcs(ii)%sid .eq. fesim%randbedingungen%mpcs(ii)%sid) then
                    fesim%internals%num_dof_vec((fesim%randbedingungen%mpcs(ii)%dependend%nid-1)*ndof+fesim%randbedingungen%mpcs(ii)%dependend%dof) = -ii
                end if
                !end do
            end if
        
        end if
    end do
end if

deallocate(fixedU)

if (textoutput .eq. .true.) write(*,*) '  ENDE  - Integriere MPCs'

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'boundary_varload_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine boundary_varload_control
