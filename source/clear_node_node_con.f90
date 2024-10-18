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
!> subroutine to delete the node list   
!
!> @details
!> this subroutine deletes the node list that is used to estimate the number
!> of non-zero elements in each row of the global stiffness matrix
!
!> @author Andreas Hauffe, TU Dresden; wissenschaftlicher Mitarbeiter 20.09.2010
!
!> $Id: clear_node_node_con.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine clear_node_node_con(virtNodes,num_nodes)
!
! Use
!
  use pre_assemble_types
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Include
!
! =================================================================================================
!
! Data types
! 
  integer, intent(in)                            :: num_nodes
  integer                                        :: ii
  
  type(virtualNode), pointer                     :: newNode    ! node which follows actNode in the list
  type(virtualNode), pointer                     :: actNode    ! active node
  type(virtualNode), dimension(num_nodes),target :: virtNodes  ! node pointer list
  
  integer                                       :: err_code=0
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
  do ii = 1, num_nodes
    actNode => virtNodes(ii)%nextNode
    DO WHILE(ASSOCIATED(actNode))
      newNode => actNode%nextNode
      if (ASSOCIATED(actNode)) deallocate(actNode)
      actNode => newNode
    END DO
  end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'clear_node_node_con'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if
  
end subroutine clear_node_node_con
