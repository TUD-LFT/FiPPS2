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
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 05.05.2010
!
!> $Id: get_node_node_con.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine get_node_node_con(fesim,virtNodes)

  use fesimulation_typen
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

!
! =================================================================================================
!
! Data types
!
  integer                                        :: ii
  integer                                        :: elem, actNodeID, nodeID
 
  type(fe_simulation)                            :: fesim
  type(virtualNode), pointer                     :: newNode, actnode
  type(virtualNode), dimension(fesim%num_nodes),target :: virtNodes
  
  integer                                        :: err_code=0
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
  ! Anlegen Vektor
  ! eine Zeile entspricht 6 Zeilen der Steifigkeitsmatrix (alle DOF eines Knotens)
  
  do ii = 1, fesim%num_nodes
    virtNodes(ii)%nodeID = ii
    virtNodes(ii)%numDofGeo = 0
    NULLIFY(virtNodes(ii)%nextNode)
  end do
  
  ! Prüfe Verbindungen zwischen Knoten von beam2-Elementen
  
  do elem = 1, size(fesim%elemente%beam2s,1)
    do actNodeID = 1, size(fesim%elemente%beam2s(elem)%nids,1)
      do nodeID = 1, size(fesim%elemente%beam2s(elem)%nids,1)
        if (nodeID .eq. actNodeID) cycle
        allocate(newNode)
        nullify(newNode%nextNode)
        newNode%nodeID = fesim%elemente%beam2s(elem)%nids(nodeID)
        newNode%numDofGeo = 6
        actNode => virtNodes(fesim%elemente%beam2s(elem)%nids(actNodeID))
        actNode%numDofGeo = 6
        if (actNode%nodeID .eq. newNode%nodeID) deallocate(newNode)
        if (ASSOCIATED(newNode)) then
          DO WHILE(ASSOCIATED(actNode%nextNode))
            actNode => actNode%nextNode
            if (actNode%nodeID .eq. newNode%nodeID) then
              deallocate(newNode)
              exit
            end if
          END DO
        end if
        if (associated(newNode)) actNode%nextNode => newNode
      end do
    end do
  end do
  
  ! Prüfe Verbindungen zwischen Knoten von quad8-Elementen
  
  do elem = 1, size(fesim%elemente%quad8s,1)
    do actNodeID = 1, size(fesim%elemente%quad8s(elem)%nids,1)
      do nodeID = 1, size(fesim%elemente%quad8s(elem)%nids,1)
        if (nodeID .eq. actNodeID) cycle
        allocate(newNode)
        nullify(newNode%nextNode)
        newNode%nodeID = fesim%elemente%quad8s(elem)%nids(nodeID)
        newNode%numDofGeo = 3
        actNode => virtNodes(fesim%elemente%quad8s(elem)%nids(actNodeID))
        actNode%numDofGeo = max(actNode%numDofGeo,3)
        if (actNode%nodeID .eq. newNode%nodeID) deallocate(newNode)
        if (ASSOCIATED(newNode)) then
          DO WHILE(ASSOCIATED(actNode%nextNode))
            actNode => actNode%nextNode
            if (actNode%nodeID .eq. newNode%nodeID) then
              deallocate(newNode)
              exit
            end if
          END DO
        end if
        if (associated(newNode)) actNode%nextNode => newNode
      end do
    end do
  end do

  ! Prüfe Verbindungen zwischen Knoten von lsolid20-Elementen
  
  do elem = 1, size(fesim%elemente%lsolid20s,1)
    do actNodeID = 1, size(fesim%elemente%lsolid20s(elem)%nids,1)
      do nodeID = 1, size(fesim%elemente%lsolid20s(elem)%nids,1)
        if (nodeID .eq. actNodeID) cycle
        allocate(newNode)
        nullify(newNode%nextNode)
        newNode%nodeID = fesim%elemente%lsolid20s(elem)%nids(nodeID)
        newNode%numDofGeo = 3
        actNode => virtNodes(fesim%elemente%lsolid20s(elem)%nids(actNodeID))
        actNode%numDofGeo = max(actNode%numDofGeo,3)
        if (actNode%nodeID .eq. newNode%nodeID) deallocate(newNode)
        if (ASSOCIATED(newNode)) then
          DO WHILE(ASSOCIATED(actNode%nextNode))
            actNode => actNode%nextNode
            if (actNode%nodeID .eq. newNode%nodeID) then
              deallocate(newNode)
              exit
            end if
          END DO
        end if
        if (associated(newNode)) actNode%nextNode => newNode
      end do
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
    write(*,*)                      'get_node_node_con'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine get_node_node_con
