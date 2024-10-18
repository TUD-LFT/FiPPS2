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
SUBROUTINE refine_mesh(verbose)

! Load modules
use konst_var
use netz_variablen

implicit none

! Declaration of variables
logical, intent(in)                               :: verbose

integer                                           :: ii, jj, i_le, iprev
integer                                           :: n_nodes_new, case_le, to_remove
double precision                                  :: x_min
double precision, dimension(:,:), allocatable     :: nodes_temp
integer, dimension(:), allocatable                :: remove_nodes

if (verbose .EQV. .TRUE.) write(*,*) 'Perform remeshing.'

! Mesh refinement using Xfoil's PANGEN routine
n_nodes_new = 300 !500 !160
allocate(nodes_temp(n_nodes_new,2))
call pangen(nodes_temp(:,1), nodes_temp(:,2), n_nodes_new, nodes(:,1), nodes(:,2), n_nodes, 0.2d0, 0.15d0, 0.2d0, 1.d0, 1.d0, 1.d0, 1.d0)
!call pangen(nodes_temp(:,1), nodes_temp(:,2), n_nodes_new, nodes(:,1), nodes(:,2), n_nodes, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, 1.d0, 1.d0)
deallocate(nodes)
allocate(nodes(n_nodes_new,2))
nodes(:,:) = nodes_temp(:,:)
deallocate(nodes_temp)
n_nodes = n_nodes_new
n_panels = n_nodes-1

x_min = minval(nodes(:,1),1)
if (x_min .lt. 0.d0) then
  nodes(:,1) = nodes(:,1) - x_min
end if

case_le = 0
i_le = minloc(nodes(:,1),1)
if ((abs(nodes(i_le,1)) .lt. 1.d-10) .and. (abs(nodes(i_le,2)) .lt. 1.d-10)) then
  case_le = 1
else
  x_min = HUGE(1.d0)
  do ii = 1, n_nodes-1
    if ((nodes(ii,2) .lt. 0.d0) .and. (nodes(ii+1,2) .gt. 0.d0) .and. (nodes(ii,1) .lt. x_min)) then
      case_le = 2
      i_le = ii
      x_min = nodes(ii,1)
    end if
  end do
end if

if (case_le .eq. 1) then
  nodes(i_le, :) = 0.d0
else if (case_le .eq. 2) then
  n_nodes_new = n_nodes + 1
  allocate(nodes_temp(n_nodes_new,2))
  nodes_temp(1:i_le,:) = nodes(1:i_le,:)
  nodes_temp(i_le+1,:) = 0.d0
  nodes_temp(i_le+2:n_nodes_new,:) = nodes(i_le+1:n_nodes,:)
  i_le = i_le + 1
  deallocate(nodes)
  allocate(nodes(n_nodes_new,2))
  nodes(:,:) = nodes_temp(:,:)
  deallocate(nodes_temp)
  n_nodes = n_nodes_new
  n_panels = n_nodes-1
else
  write(*,*) 'ERROR in adding point (0,0).'
  STOP
end if

to_remove = 0
allocate(remove_nodes(n_nodes))
remove_nodes(:) = -1
iprev = 1
do ii = 2, i_le
  if (nodes(ii,1) .ge. nodes(iprev,1)) then
    to_remove = to_remove + 1
    remove_nodes(to_remove) = ii
  else
    iprev = ii
  end if
end do

iprev = i_le
do ii = i_le+1, n_nodes
  if (nodes(ii,1) .le. nodes(iprev,1)) then
    to_remove = to_remove + 1
    remove_nodes(to_remove) = ii
  else
    iprev = ii
  end if
end do

if (to_remove .gt. 1.d-2*n_nodes) then
  write(*,*) '  Error performing remesh!'
  write(*,*) '  More than 1% of the nodes have to be removed.'
  stop
end if

if (to_remove .gt. 0) then
  n_nodes_new = n_nodes - to_remove
  allocate(nodes_temp(n_nodes_new,2))
  jj = 1
  do ii = 1, n_nodes
    if (ii .ne. remove_nodes(jj)) then
      nodes_temp(ii-jj+1,:) = nodes(ii,:)
    else
      jj = jj+1
    end if
  end do
  deallocate(nodes)
  allocate(nodes(n_nodes_new,2))
  nodes(:,:) = nodes_temp(:,:)
  deallocate(nodes_temp)
  n_nodes = n_nodes_new
  n_panels = n_nodes-1
end if  
deallocate(remove_nodes)

END SUBROUTINE refine_mesh
