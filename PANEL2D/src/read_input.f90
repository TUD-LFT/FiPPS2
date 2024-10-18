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
SUBROUTINE read_input(filename,verbose)

! Load modules
use konst_var
use netz_variablen
use functions

implicit none

! Declaration of variables
intent(in)                                        :: filename,verbose
character(len=*)                                  :: filename
character(len=256)                                :: af_file

logical                                           :: delta_aoa
logical                                           :: verbose

double precision                                  :: aoa_diff
integer                                           :: ii, jj
integer                                           :: unit_in=21
integer                                           :: io_error

! Calculation
if (verbose) write(*,*) '  BEGIN: Read input data.'

! Open input file
open(unit=unit_in,file=TRIM(filename),status='old',action='read',iostat=io_error)
if (io_error .EQ. 0) then
    read(unit_in,*)
    read(unit_in,'(L)') cl_given
    read(unit_in,*)
    read(unit_in,'(L)') delta_aoa
    read(unit_in,*)
    if (delta_aoa .EQV. .FALSE.) then
        read(unit_in,*) aoa
    else
        read(unit_in,*) aoa_diff
        aoa = rad_to_deg(aoa) + aoa_diff
    end if
    read(unit_in,*)
    read(unit_in,*) cl_target
    read(unit_in,*)
    read(unit_in,*) reynolds
    read(unit_in,*)
    read(unit_in,*) Ma
    read(unit_in,*)
    read(unit_in,*) q_infty
    read(unit_in,*)
    read(unit_in,*) chordscale
    read(unit_in,*)
    read(unit_in,*) af_file
    read(unit_in,*)
    read(unit_in,*) res_file
else
    write(*,*)
    write(*,*) '  ERROR: File ', TRIM(filename), ' cannot be read!'
    write(*,*)
end if
! Close input file
close(unit_in)

! Open file with airfoil coordinates
open(unit=unit_in,file=TRIM(af_file),status='old',action='read',iostat=io_error)
if (io_error .EQ. 0) then
    read(unit_in,*)
    read(unit_in,*) n_nodes
    allocate(nodes(n_nodes,2))
    read(unit_in,*)
    do ii = 1, n_nodes
        read(unit_in,*) (nodes(ii,jj), jj=1,2)
    end do
    n_panels = n_nodes-1
else
    write(*,*)
    write(*,*) '  ERROR: File ', TRIM(af_file), ' cannot be read!'
    write(*,*)
end if
! Close file with airfoil coordinates
close(unit_in)

! Transform angle of attack to radiant
aoa = deg_to_rad(aoa)

! Sort nodes in inverse order, if necessary
call sort_nodes()

if (verbose) write(*,*) '  END:   Read input data.'

END SUBROUTINE read_input
