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
! Copyright (C) 2010-2012  Daniel Filkovic

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %                                                      %
! %            APAME - Aircraft Panel Method             %
! %______________________________________________________%
! %                                                      %
! %              3D potential flow solver                %
! %                                                      %
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This file is part of APAME.

! APAME is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! APAME is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with APAME.  If not, see <http://www.gnu.org/licenses/>.

! file apame.f90

program apame

implicit none

interface
    subroutine apame_main(fname,readmodel,displacements,dispfac,pressures,F_total,cdi,output)
        use module_kind_and_konst, only: kind_float
        character(len=*), optional                                                 :: fname             !< Filename of the apame input file
        logical,                                            intent(in), optional   :: readmodel         !< Flag, to read the input file
        real(kind=kind_float), dimension(:,:),              intent(in), optional   :: displacements     !< Displacements to be added on node coordinates
        real(kind=kind_float),                              intent(in), optional   :: dispfac           !< Factor with which displacements are added
        real(kind=kind_float), dimension(:,:), allocatable, intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        real(kind=kind_float), dimension(:), allocatable,   intent(out), optional  :: F_total
        real(kind=kind_float), dimension(:), allocatable,   intent(out), optional  :: cdi
        logical,                                            intent(in), optional   :: output            !< Flag, if APAME should write any output
    end subroutine apame_main
end interface

call apame_main()

end program apame
