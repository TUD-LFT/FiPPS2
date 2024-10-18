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
!
SUBROUTINE vtkoutput(Utot)

  USE globale_variablen

  IMPLICIT NONE
  INTEGER :: ii, jj, kk, ll
  INTEGER :: cell_number, cell_data_number
  
  integer, dimension(3)	:: tria3_node_ids
  integer, dimension(4)	:: quad4_node_ids
  
  double precision, dimension(num_dof,num_subcases), intent(in) :: Utot

  OPEN(UNIT=42,FILE='output.vtk')

  WRITE(42, '(A26)') '# vtk DataFile Version 4.0'
  WRITE(42, '(A12)') 'Cube example'
  WRITE(42, '(A5)')  'ASCII'
  WRITE(42, '(A25)') 'DATASET UNSTRUCTURED_GRID'
  WRITE(42, '(A7,I8,A7)') 'POINTS ', size(nodes,1), ' double'
  
  DO ii = 1, size(nodes,1)
      WRITE(42, '(E24.17,X,E24.17,X,E24.17)') nodes(ii)%coords(1), nodes(ii)%coords(2), nodes(ii)%coords(3)
  END DO
  
  cell_number      = 0
  cell_data_number = 0
  
  IF (ALLOCATED(tria3s) .EQ. .TRUE.) THEN 
    cell_number      = cell_number + size(tria3s,1)
    cell_data_number = cell_data_number + size(tria3s,1) * 4
  END IF
  
  IF (ALLOCATED(quad4s) .EQ. .TRUE.) THEN 
    cell_number      = cell_number + size(quad4s,1)
    cell_data_number = cell_data_number + size(quad4s,1) * 5
  END IF
  
  WRITE(42, '(A6,I8,X,I8)') 'CELLS ', cell_number, cell_data_number
  
  IF (ALLOCATED(tria3s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(tria3s,1)
    
      ! get external node ids
      do kk=1,size(nodes,1)
        do ll=1,3
          if (tria3s(ii)%nids(ll) == nodes(kk)%nid) then
      
            ! get internal node ids
            tria3_node_ids(ll) = kk
         
          end if
        end do
      end do
      
      WRITE(42, '(3(I10,X),I10)') 3, tria3_node_ids(1)-1, tria3_node_ids(2)-1, tria3_node_ids(3)-1
    END DO
  END IF
  
  IF (ALLOCATED(quad4s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(quad4s,1)
    
      ! get external node ids
      do kk=1,size(nodes,1)
        do ll=1,4
          if (quad4s(ii)%nids(ll) == nodes(kk)%nid) then
      
            ! get internal node ids
            quad4_node_ids(ll) = kk
         
          end if
        end do
      end do
      
      WRITE(42, '(4(I10,X),I10)') 4, quad4_node_ids(1)-1, quad4_node_ids(2)-1, quad4_node_ids(3)-1, quad4_node_ids(4)-1
    END DO
  END IF
 
  WRITE(42, '(A11,I8)') 'CELL_TYPES ', cell_number
  
  IF (ALLOCATED(tria3s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(tria3s,1)
      WRITE(42, '(I10)') 5
    END DO
  END IF
  
  IF (ALLOCATED(quad4s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(quad4s,1)
      WRITE(42, '(I10)') 9
    END DO
  END IF

  WRITE(42, '(A11,I8)') 'POINT_DATA ', size(nodes,1)
  DO kk = 1,num_subcases

    WRITE(42, '(A17,I8.8,A7)') 'VECTORS Lastfall_', kk,' double'
    
    DO jj = 1, size(nodes,1)
        WRITE(42, '(E17.10,X,E17.10,X,E17.10)') Utot((jj-1)*6+1,kk), Utot((jj-1)*6+2,kk), Utot((jj-1)*6+3,kk)
    END DO
    
  END DO
  
  CLOSE(42)

END SUBROUTINE vtkoutput
