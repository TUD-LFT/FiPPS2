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
PROGRAM PANEL2D

implicit none

character(len=256)                                                                 :: input_file
double precision, dimension(:), allocatable                                        :: Mcrit
 
interface
    subroutine panel2d_main(input_file, readmodel, inviscid, displacements, dispfac, given_dalpha, paneling, pressures, F_total, al, ca, cd, cm, hkmax, xsep, Mcrit, output)
        character(len=*),                                   intent(in)             :: input_file        !< Filename of the apame input file
        logical,                                            intent(in),  optional  :: readmodel         !< Flag, to read the input file
        logical,                                            intent(in),  optional  :: inviscid          !< Flag, to force inviscid calculation
        double precision, dimension(:,:),                   intent(in),  optional  :: displacements     !< Displacements to be added on node coordinates
        double precision,                                   intent(in),  optional  :: dispfac           !< Factor with which displacements are added
        double precision,                                   intent(in),  optional  :: given_dalpha      !< Delta Alpha (overrides value from inputfile)
        logical,                                            intent(in),  optional  :: paneling          !< Enable repaneling in PANEL2D
        double precision, dimension(:,:), allocatable,      intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        double precision, dimension(:), allocatable,        intent(out), optional  :: F_total
        double precision, dimension(:), allocatable,        intent(out), optional  :: al
        double precision, dimension(:), allocatable,        intent(out), optional  :: ca
        double precision, dimension(:), allocatable,        intent(out), optional  :: cd
        double precision, dimension(:), allocatable,        intent(out), optional  :: cm
        double precision, dimension(:), allocatable,        intent(out), optional  :: hkmax
        double precision, dimension(:), allocatable,        intent(out), optional  :: xsep
        double precision, dimension(:), allocatable,        intent(out), optional  :: Mcrit
        logical,                                            intent(in),  optional  :: output            !< Flag, if PANEL2D should write any output
    end subroutine panel2d_main
end interface

! Inputfiles
call get_command_argument(1,input_file)
if (trim(input_file) .eq. '') then
    write(*,*) 'No input-file specified.'
    STOP
end if

call panel2d_main(input_file,output=.TRUE.,paneling=.TRUE.,Mcrit=Mcrit)

END PROGRAM PANEL2D
