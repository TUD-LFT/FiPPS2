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
SUBROUTINE xfw_read_input(fname)

use xfw_konst_var
use xfw_netz_variablen
use xfw_functions

implicit none

! Variablendeklaration
intent(in)                                        :: fname
character(len=*)                                  :: fname
character(len=256)                                :: af_file

double precision                                  :: z_side1, z_side2, aoa_diff
double precision, dimension(:,:), allocatable     :: nodes_temp
integer                                           :: le_pos
integer                                           :: ii, jj
integer                                           :: unit_in=21
integer                                           :: io_error

! Berechnung

! Oeffne Eingabedatei
open(unit=unit_in,file=TRIM(fname),status='old',action='read',iostat=io_error)
if (io_error .EQ. 0) then
    read(unit_in,*)
    read(unit_in,'(L)') cl_given
    read(unit_in,*)
    read(unit_in,'(L)') delta_aoa
    read(unit_in,*)
    if (delta_aoa .EQV. .FALSE.) then
        read(unit_in,'(E24.16)') aoa
    else
        read(unit_in,'(E24.16)') aoa_diff
        aoa = aoa + aoa_diff
    end if
    read(unit_in,*)
    read(unit_in,'(E24.16)') cl_target
    read(unit_in,*)
    read(unit_in,'(E24.16)') reynolds
    read(unit_in,*)
    read(unit_in,'(E24.16)') Ma
    read(unit_in,*)
    read(unit_in,'(E24.16)') q_infty
    read(unit_in,*)
    read(unit_in,'(E24.16)') chordscale
    read(unit_in,*)
    read(unit_in,*) af_file
    read(unit_in,*)
    read(unit_in,*) res_file
    read(unit_in,*)
    read(unit_in,'(A256)') start_xfoil
else
    write(*,*)
    write(*,*) '  FEHLER: Datei ', TRIM(fname), ' kann nicht gelesen werden!'
    write(*,*)
end if
! Schliesse Eingabedatei
close(unit_in)

! Oeffne Datei mit Profilkoordinaten
open(unit=unit_in,file=TRIM(af_file),status='old',action='read',iostat=io_error)
if (io_error .EQ. 0) then
    read(unit_in,*)
    read(unit_in,'(I10)') n_nodes
    allocate(nodes(n_nodes,2))
    read(unit_in,*)
    do ii = 1, n_nodes
        read(unit_in,'(2E24.16)') (nodes(ii,jj), jj=1,2)
    end do
    n_panels = n_nodes-1
else
    write(*,*)
    write(*,*) '  FEHLER: Datei ', TRIM(af_file), ' kann nicht gelesen werden!'
    write(*,*)
end if
! Schliesse Datei mit Profilkoordinaten
close(unit_in)

! Bei Bedarf Sortierung der Knoten umdrehen
! Fuer die Berechnung muessen die Knoten von der Hinterkante ueber
! die Oberseite, die Unterseite zurueck zur Hinterkante laufen
! Finde Knotenposition der Vorderkante
le_pos  = minloc(nodes(:,1), 1)
! Pruefe, ob Knotennummerierung umgedreht werden muss
z_side1 = sum(nodes(1:le_pos,2))
z_side2 = sum(nodes(le_pos+1:,2))
if (z_side2 .gt. z_side1) then
    invert_points = .true.
else
    invert_points = .false.
end if

if (invert_points .eqv. .true.) then
    ! Drehe Knotensortierung um
    allocate(nodes_temp(n_nodes,2))
    nodes_temp(:,:) = nodes(:,:)
    do ii = 1, n_nodes
        nodes(ii,:) = nodes_temp(n_nodes-ii+1,:)
    end do
    deallocate(nodes_temp)
end if

END SUBROUTINE xfw_read_input
