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
!> @brief
!> Module
!
!> @details
!
!> @author
!
!
!> $Id: dof_types.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module dof_types


  implicit none
  
  integer, parameter :: type_free_dof          = 1      ! freien Freiheiten
  integer, parameter :: type_constrained_dof   = 2      ! vorgegebenen Freiheiten (inhomogene Randbedingungen)
  integer, parameter :: type_slave_coupled_dof = 3      ! Slave-Freiheiten von Couplings
  integer, parameter :: type_fixed_dof         = 4      ! gesperrten Freiheiten

  integer            :: num_free_dof                    ! Anzahl der freien Freiheiten
  integer            :: num_constrained_dof             ! Anzahl der vorgegebenen Freiheiten (inhomogene Randbedingungen)
  integer            :: num_slave_coupled_dof           ! Anzahl der Slave-Freiheiten von Couplings
  integer            :: num_fixed_dof                   ! Anzahl der gesperrten Freiheiten
  
  integer, dimension(1:4) :: start_type_dof             ! 

end module dof_types
