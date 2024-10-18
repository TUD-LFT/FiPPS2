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
SUBROUTINE write_vtk(vtk_filename,scaleFac,cp,ue,uen,p_mano,bl_values)

use konst_var
use netz_variablen

implicit none

! Variablendeklaration
intent(in)                                        :: vtk_filename,scaleFac,cp,ue,uen,p_mano,bl_values

character(len=*)                                  :: vtk_filename
double precision                                  :: scaleFac
double precision, dimension(2)                    :: vec_temp
double precision, dimension(n_panels)             :: cp, ue, p_mano
double precision, dimension(n_nodes)              :: uen
type(bl_point), dimension(n_nodes)                :: bl_values

integer                                           :: ii,jj
integer                                           :: unit_id = 31

open(unit=unit_id,file=TRIM(vtk_filename),status='replace')

! Write header
write(unit_id,'(A26)') '# vtk DataFile Version 4.0'
write(unit_id,'(A7)')  'PANEL2D'
write(unit_id,'(A5)')  'ASCII'
write(unit_id,'(A25)') 'DATASET UNSTRUCTURED_GRID'

! Write nodes
write(unit_id,'(A6,X,I0,X,A6)') 'POINTS', n_nodes, "DOUBLE"
do ii = 1, n_nodes
    write(unit_id,'(3E24.16)') (nodes(ii,jj)*scaleFac, jj=1,2), 0.d0
end do

! Write panels
write(unit_id,'(A5,X,I0,X,I0)') 'CELLS', n_panels, 3*n_panels
do ii = 1, n_panels
    write(unit_id,'(I0,X,I0,X,I0)') 2, ii-1, ii
end do
write(unit_id,'(A10,X,I0)') 'CELL_TYPES', n_panels
do ii = 1, n_panels
    write(unit_id,'(I0)') 3
end do

! Write cell data
write(unit_id,'(A9,X,I0)') 'CELL_DATA', n_panels
write(unit_id,'(A7,X,A3,X,A6)') 'SCALARS', 'C_P', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_panels
    write(unit_id,'(E24.16)') cp(ii)
end do
write(unit_id,'(A7,X,A6,X,A6)') 'SCALARS', 'P_MANO', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_panels
    write(unit_id,'(E24.16)') p_mano(ii)
end do
write(unit_id,'(A7,X,A2,X,A6)') 'SCALARS', 'UE', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_panels
    write(unit_id,'(E24.16)') ue(ii)
end do
write(unit_id,'(A7,X,A6,X,A6)') 'NORMALS', 'NORMAL', 'DOUBLE'
do ii = 1, n_panels
    write(unit_id,'(3E24.16)') (n_vec(ii,jj), jj=1,2), 0.d0
end do
write(unit_id,'(A7,X,A5,X,A6)') 'NORMALS', 'T_VEC', 'DOUBLE'
do ii = 1, n_panels
    write(unit_id,'(3E24.16)') (t_vec(ii,jj), jj=1,2), 0.d0
end do

! Write point data
write(unit_id,'(A10,X,I0)') 'POINT_DATA', n_nodes
write(unit_id,'(A7,X,A3,X,A6)') 'SCALARS', 'H_K', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_nodes
    write(unit_id,'(E24.16)') bl_values(ii)%H
end do
write(unit_id,'(A7,X,A3,X,A6)') 'SCALARS', 'C_F', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_nodes
    write(unit_id,'(E24.16)') bl_values(ii)%cf
end do
write(unit_id,'(A7,X,A6,X,A6)') 'SCALARS', 'DELTAS', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_nodes
    write(unit_id,'(E24.16)') bl_values(ii)%deltas
end do
write(unit_id,'(A7,X,A5,X,A6)') 'SCALARS', 'THETA', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_nodes
    write(unit_id,'(E24.16)') bl_values(ii)%theta
end do
write(unit_id,'(A7,X,A4,X,A6)') 'SCALARS', 'NFAC', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_nodes
    write(unit_id,'(E24.16)') bl_values(ii)%nfac
end do
write(unit_id,'(A7,X,A2,X,A6)') 'SCALARS', 'UE', 'DOUBLE'
write(unit_id,'(A12,X,A7)') 'LOOKUP_TABLE', 'default'
do ii = 1, n_nodes
    write(unit_id,'(E24.16)') uen(ii)
end do
write(unit_id,'(A7,X,A12,X,A6)') 'VECTORS', 'BL_THICKNESS', 'DOUBLE'
do ii = 1, n_nodes
    if (ii .eq. 1) then
        vec_temp(:) = n_vec(1,:)
    else if (ii .eq. n_nodes) then
        vec_temp(:) = n_vec(n_panels,:)
    else
        vec_temp(:) = (n_vec(ii,:) + n_vec(ii-1,:))/2.d0
    end if
    vec_temp(:) = vec_temp(:)*bl_values(ii)%deltas
    write(unit_id,'(3E24.16)') (vec_temp(jj), jj=1,2), 0.d0
end do

! close file
close(unit_id)

END SUBROUTINE write_vtk
