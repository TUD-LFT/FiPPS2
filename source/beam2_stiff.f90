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
!> subroutine for calculation of 2-node beam element stiffness matrix
!
!> @details 
!> subroutine calculates the stiffness matrix of a 2-node straight beam element
!> with constant cross-section properties
!
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 22.06.2010
!
!> $Id: beam2_stiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_stiff(fesim,ii,beam2_node_coords,AA,I11,I22,It,Emat,Gmat,beam2_K,le)
!
! Use
!
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
type(fe_simulation), intent(in)                     :: fesim
integer, intent(in)                                 :: ii                !< beam element index
double precision, dimension(2,3), intent(in)        :: beam2_node_coords !< nodal coordinates of element nodes
double precision, intent(in)                        :: AA    !< beam cross section area
double precision, intent(in)                        :: I11   !< beam cross section moments of inertia  about z = Izz
double precision, intent(in)                        :: I22   !< beam cross section moments of inertia about y = Iyy
double precision, intent(in)                        :: It    !< beam cross section torsional moment of inertia
double precision, intent(in)                        :: Emat  !< material youngs
double precision, intent(in)                        :: Gmat  !< material shear modulus
!
! Output
!
double precision, intent(out)                       :: le       !< element length
!
! Input + Output
!
double precision, dimension(12,12), intent(inout)   :: beam2_K  !< beam2-element stiffness matrix
!
! Inner
!
double precision                                    :: facDP
double precision, dimension(3)                      :: vv
double precision, dimension(12,12)                  :: TRMat
integer                                             :: jj,kk
integer                                             :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
facDP = 0.0D0
beam2_K = 0.0D0
!
! =================================================================================================
!
! Calculation

  ! get vector v (see MSC.NASTRAN 2005 Reference Guide p. 59)
  
  !vv = (beam2s(ii)%xi(jj), jj=1,3)
  vv(1:3)=fesim%elemente%beam2s(ii)%xi(1:3)

  ! get transformation matrix from local to global coordinate system
  
  ! call rotation_beam2 (beam2_node_coords,(beam2s(ii)%xi(jj), jj=1,3),'lg',TRMat,le)
   call beam2_rotation (beam2_node_coords,vv,'lg',TRMat,le)
  
  ! fill upper diagonal of element stiffness matrix
  
  beam2_K(1,1)   =  Emat*AA/le                  ! 1st row
  beam2_K(1,7)   = -beam2_K(1,1)
  beam2_K(2,2)   =  12.D0*Emat*I11/(le**3)      ! 2nd row
  beam2_K(2,6)   =  6.D0*Emat*I11/(le**2)
  beam2_K(2,8)   = -beam2_K(2,2)
  beam2_K(2,12)  =  beam2_K(2,6)
  beam2_K(3,3)   =  12.D0*Emat*I22/(le**3)      ! 3rd row
  beam2_K(3,5)   = -6.D0*Emat*I22/(le**2)
  beam2_K(3,9)   = -beam2_K(3,3)
  beam2_K(3,11)  =  beam2_K(3,5)
  beam2_K(4,4)   =  Gmat*It/le                  ! 4th row
  beam2_K(4,10)  = -beam2_K(4,4)
  beam2_K(5,5)   =  4.D0*Emat*I22/le            ! 5th row
  beam2_K(5,9)   =  beam2_K(5,5)*1.5D0/le
  beam2_K(5,11)  =  beam2_K(5,5)/2.0D0
  beam2_K(6,6)   =  4.D0*Emat*I11/le            ! 6th row
  beam2_K(6,8)   = -beam2_K(6,6)*1.5D0/le
  beam2_K(6,12)  =  beam2_K(6,6)/2.0D0
  beam2_K(7,7)   =  beam2_K(1,1)                ! 7th row
  beam2_K(8,8)   =  beam2_K(2,2)                ! 8th row
  beam2_K(8,12)  =  beam2_K(6,8)
  beam2_K(9,9)   =  beam2_K(3,3)                ! 9th row
  beam2_K(9,11)  =  beam2_K(5,9)
  beam2_K(10,10) =  beam2_K(4,4)                ! 10th row
  beam2_K(11,11) =  beam2_K(5,5)                ! 11th row
  beam2_K(12,12) =  beam2_K(6,6)                ! 12th row
  
  ! fill lower diagonal
  
  do jj=2,12
    do kk=1,(jj-1)
      beam2_K(jj,kk) = beam2_K(kk,jj)
    end do
  end do
  
  ! transform element stiffness matrix to global coordinates
  
  beam2_K = matmul(matmul(TRMat,beam2_K),transpose(TRMat))
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_stiff'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine beam2_stiff
