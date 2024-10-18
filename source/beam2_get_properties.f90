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
!> Subroutine emits the section and isotropic material properties for a 2-node
!> beam element
!
!> @details
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 24.06.2010
!
!> $Id: beam2_get_properties.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_get_properties(fesim,ii,AA,I11,I22,I12,It,Emat,Gmat)

  use fesimulation_typen
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
! Input
!
integer, intent(in)                             :: ii    !< Beam2-Element-ID
!
! Output
!
double precision, intent(out)                   :: AA    !< beam cross section area
double precision, intent(out)                   :: I11   !< beam cross section moments of inertia about
double precision, intent(out)                   :: I22   !< beam cross section moments of inertia about
double precision, intent(out)                   :: I12   !< beam cross section product of inertia
double precision, intent(out)                   :: It    !< beam cross section torsional moment of inertia
double precision, intent(out)                   :: Emat  !< material youngs modulus
double precision, intent(out)                   :: Gmat  !< material shear modulus
!
! Input + Output
!
type(fe_simulation), intent(inout)              :: fesim
!
! Internal
!
double precision                                :: nu
double precision                                :: theta,Ixi, Ieta, Ietaxi

integer                                         :: kk
integer                                         :: err_code=0
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
  ! set beam section properties
  
  AA    = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%AA
  I11   = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%I11
  I22   = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%I22
  I12   = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%I12
  
  ! get beam material stiffness
  
  do kk=1,size(fesim%materialien%mat1s,1)
     if (fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%mid == fesim%materialien%mat1s(kk)%mid) then
        
        call mat1_calc_missing_value(fesim,kk,Emat,Gmat,nu)
        fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%intMat1ID = kk
        
     end if
  end do
  
  ! angle between principal axes and element local coordinate system
  
  if (fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%angle /= 0.D0) then

    theta=fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%angle

    Ieta = (I22+I11)/2.D0+(I22-I11)/2.D0*cos(2.D0*theta)*-I12*sin(2.D0*theta)
    Ixi  = (I22+I11)/2.D0-(I22-I11)/2.D0*cos(2.D0*theta)*-I12*sin(2.D0*theta)

    Ietaxi = (I22-I11)/2.D0*sin(2.D0*theta)*+I12*cos(2.D0*theta)

    I11 = Ixi
    I22 = Ieta
    I12 = Ietaxi

  end if
  
  ! determine torsional moment of inertia
  
  if (fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%It /= 0.D0) then
    It  = fesim%eigenschaften%pbeams(fesim%elemente%beam2s(ii)%int_pid)%It
  else
    ! if no torsional moment of inertia is specified, the polar moment of inertia
    ! is used instead; this is done according ANSYS (ANSYS Help -> Element Reference
    ! -> Element Library -> BEAM4 -> BEAM4 Input Data)
    ! the torsional moment of inertia is usually less than the polar moment of inertia
    !
    ! for solution type 600 (sol600) in NASTRAN It is taken as 0.5*(I11+I22)
    ! (NASTRAN Quick Reference Guide -> Bulk Data Entries -> PBEAM)

    It  = I11+I22                        ! polar moment of inertia
  end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_get_properties'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine beam2_get_properties
