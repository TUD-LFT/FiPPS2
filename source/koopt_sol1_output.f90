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
!> $Id: koopt_sol1_output.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine koopt_sol1_output(stress, shell_strain_top, shell_strain_bot, Utot, pbeamsize, &
                             umax, beam_stresses, shell_strains_max, shell_strains_min, &
                             num_pshells, num_parts)

use globale_variablen
use netz_variablen
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
! Input
!
  double precision, dimension(num_elements, 7), intent(in) :: stress 
  double precision, dimension(num_elements, 3), intent(in) :: shell_strain_top, shell_strain_bot
  double precision, dimension(num_dof), intent(in)         :: Utot
  integer                                                  :: pbeamsize
  integer                                                  :: num_parts, num_pshells
!
! Internal
!
  integer                                                  :: elem, node, ii, offset
  double precision                                         :: absvec, umax
  double precision, dimension(pbeamsize,4)                 :: beam_stresses
  double precision, dimension(num_parts,3)                 :: shell_strains_max,shell_strains_min

  integer                                                  :: err_code=0
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
  beam_stresses(:,1) = -HUGE(1.d0)
  beam_stresses(:,2) =  HUGE(1.d0)
  beam_stresses(:,3) = -HUGE(1.d0)
  beam_stresses(:,4) = -HUGE(1.d0)
  
  do elem = 1, size(beam2s,1)
  
    ! maximale Spannung
    if (stress(beam2s(elem)%eid,4) .GT. beam_stresses(beam2s(elem)%int_pid,1)) then
      beam_stresses(beam2s(elem)%int_pid,1) = stress(beam2s(elem)%eid,4)
    end if
    
    ! minimale Spannung
    if (stress(beam2s(elem)%eid,5) .LT. beam_stresses(beam2s(elem)%int_pid,2)) then
      beam_stresses(beam2s(elem)%int_pid,2) = stress(beam2s(elem)%eid,5)
    end if
    
    ! maximale L�ngskraft
    if (stress(beam2s(elem)%eid,6) .GT. beam_stresses(beam2s(elem)%int_pid,3)) then
      beam_stresses(beam2s(elem)%int_pid,3) = stress(beam2s(elem)%eid,6)
    end if
    
    ! maximales Biegemoment
    if (stress(beam2s(elem)%eid,7) .GT. beam_stresses(beam2s(elem)%int_pid,4)) then
      beam_stresses(beam2s(elem)%int_pid,4) = stress(beam2s(elem)%eid,7)
    end if
  
  end do
  
  umax = 0.d0
  do node = 1, num_nodes
    absvec = Utot((node-1)*ndof+1) * Utot((node-1)*ndof+1) + &
             Utot((node-1)*ndof+2) * Utot((node-1)*ndof+2) + &
             Utot((node-1)*ndof+3) * Utot((node-1)*ndof+3)
    absvec = sqrt(absvec)
    if (absvec .GT. umax) umax = absvec
  end do
  
  shell_strains_max(:,1:3) = -1.d300
  shell_strains_min(:,1:3) =  1.d300
  
  do elem = 1, size(quad8s,1)
   
    if (quad8s(elem)%propType == 1) then
      offset = 0
    else if (quad8s(elem)%propType == 2) then
      offset = num_pshells
    end if

    do ii = 1,3
      if (shell_strain_top(quad8s(elem)%eid,ii) .GT. shell_strains_max(quad8s(elem)%int_pid+offset,ii)) then
        shell_strains_max(quad8s(elem)%int_pid+offset,ii) = shell_strain_top(quad8s(elem)%eid,ii)
      end if
      if (shell_strain_bot(quad8s(elem)%eid,ii) .GT. shell_strains_max(quad8s(elem)%int_pid+offset,ii)) then
        shell_strains_max(quad8s(elem)%int_pid+offset,ii) = shell_strain_bot(quad8s(elem)%eid,ii)
      end if
      
      if (shell_strain_top(quad8s(elem)%eid,ii) .LT. shell_strains_min(quad8s(elem)%int_pid+offset,ii)) then
        shell_strains_min(quad8s(elem)%int_pid+offset,ii) = shell_strain_top(quad8s(elem)%eid,ii)
      end if
      if (shell_strain_bot(quad8s(elem)%eid,ii) .LT. shell_strains_min(quad8s(elem)%int_pid+offset,ii)) then
        shell_strains_min(quad8s(elem)%int_pid+offset,ii) = shell_strain_bot(quad8s(elem)%eid,ii)
      end if
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
   write(*,*)                      'koopt_sol1_output'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine koopt_sol1_output
