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
!> @author 
!
!> $Id: sol2_output.f90 414 2019-10-17 11:45:47Z DOM+ahauffe $
!> $Author: DOM+ahauffe $
!> $Revision: 414 $
!> $Date: 2019-10-17 13:45:47 +0200 (Do, 17. Okt 2019) $
!
! =================================================================================================
SUBROUTINE sol1_output_control(fesim,Utot,scloop,cdi,aeroConverged)
!
! Use
!
  use fesimulation_typen
!
! =================================================================================================
!
  IMPLICIT NONE
  
  INTERFACE 
    SUBROUTINE SOL1_OUTPUT_GLAWI2(FESIM,UTOT,SCLOOP,CDI,AEROCONVERGED)
      USE FESIMULATION_TYPEN
      TYPE (FE_SIMULATION) :: FESIM
      REAL(KIND=8) :: UTOT(FESIM%NUM_DOF)
      INTEGER(KIND=4), INTENT(IN) :: SCLOOP
      REAL(KIND=8) :: CDI(:)
      LOGICAL(KIND=4) :: AEROCONVERGED
    END SUBROUTINE SOL1_OUTPUT_GLAWI2
  END INTERFACE 
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
! Input
!
  type(fe_simulation), intent(inout)                            :: fesim
  double precision, dimension(fesim%num_dof), intent(in)        :: Utot
  integer, intent(in)                                           :: scloop
  double precision, dimension(:), allocatable, intent(in)       :: cdi
  logical, intent(in)                                           :: aeroConverged
!
! Internal
!
  integer                                                       :: err_code=0
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
  !
  !     #ifdef koopt_output
  !         call sol1_output_koopt(Utot)
  !     #endif
  if (fesim%ausgabe%outputVTK .eq. .true.) then
      call sol1_output_vtk(fesim,Utot,scloop)
  end if
  if (fesim%ausgabe%outputShort .eq. .true.) then
      call sol1_output_sab(fesim,Utot,scloop)
  endif
  !     #ifdef optitube_output
  !         call sol1_output_optitube(Utot)
  !     #endif
  if (fesim%ausgabe%outputGlawi .eq. .true.) then
      call sol1_output_glawi2(fesim,Utot,scloop,cdi,aeroConverged)
  end if
  if (fesim%ausgabe%outputOptitube .eq. .true.) then
      call sol1_output_glawi(fesim,Utot,scloop)
  end if
  if (fesim%ausgabe%outputHyMoWi .eq. .true.) then
      call sol1_output_hymowi2(fesim,Utot,scloop,aeroConverged)
  end if
  if (fesim%ausgabe%outputAdviLa .eq. .true.) then
      call sol1_output_advila(fesim,Utot,scloop)
  end if
  !     #ifdef hymowi_output
  !         call sol1_output_hymowi(Utot)
  !     #endif
  !
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                       'sol1_output_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE sol1_output_control
