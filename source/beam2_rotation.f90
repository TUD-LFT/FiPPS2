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
!> Subroutine for transformation matrix for beam2-element
!
!> @details 
!> Subroutines computes the transformation matrix from local to global
!> coordinates (or the other way around) for a 2-node beam element
!> Additionally he element length is computed
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 22.06.2010
!
!> $Id: beam2_rotation.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
! 
! =================================================================================================
subroutine beam2_rotation (beam2_node_coords,vv,mode,TRMat,le)

use konstanten
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Interface
!

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
double precision, dimension (2,3), intent(in)    :: beam2_node_coords   !< nodal coordinates of element nodes in global coordinates
double precision, dimension(3), intent(in)       :: vv                  !< vector in x-y-plane of local element coordinate system                          
character(2), intent(in)                         :: mode                !< determines whether transformation matrix is for transormation local to global or global to local system
!
! Output
!
double precision, dimension (12,12), intent(out) :: TRMat               !< transformation matrix
double precision, intent(out)                    :: le                  !< element length
!
! Input + Output
!

!
! Internal
!
double precision, dimension(3,3)                 :: BC
integer                                          :: kk,ll,mm
integer                                          :: imat,jmat

integer                                          :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
  TRmat = 0.d0
!
! =================================================================================================
!
! Calculation
!
  !
  ! calculation of cosine matrix
  !
  
  call beam2_cosinematrix(beam2_node_coords,vv,BC,le)

  !
  ! blow up to transformation matrix
  !

  do kk=1,(2*ndof/3)    ! 2*3 Translations, 2*3 Rotations
    
    imat=(kk-1)*3
    jmat=(kk-1)*3

    do ll=1,3
      do mm=1,3
        
        TRMat(imat+ll,jmat+mm) = BC(ll,mm)
        
      end do
    end do
    
  end do

  ! transformation matrix so far is for trafo local -> global

  if (mode == 'lg') then
  
  ! do nothing

  else if (mode == 'gl') then
  
    TRMat = transpose(TRMat)
   
  else
   
    write(*,*) 'Error computing transformation matrix'
    write(*,*) 'Input .mode. falsely defined'
    err_code = 1
    goto 9999
   
  end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_rotation'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine beam2_rotation
