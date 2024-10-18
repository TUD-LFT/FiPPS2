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
!> subroutine for calculation of 2-node beam element geometric stiffness matrix
!
!> @details 
!> subroutine calculates the stiffness matrix of a 2-node straight beam element
!> with constant cross-section properties
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter, 22.06.2010
!> @author Florian Dexl, TU Dresden, wissenschaftlicher Mitarbeiter, 12.04.2019
!
!> $Id: beam2_geostiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine beam2_geostiff(fesim,ii,disp_beam2,beam2_node_coords,AA,I11,I22,Emat,btemp,alpha,beam2_Kg)
!
use fesimulation_typen
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
type(fe_simulation)                                 :: fesim
integer, intent(in)                                 :: ii                   !< beam element index
double precision, dimension(2,3), intent(in)        :: beam2_node_coords    !< nodal coordinates of element nodes
double precision, intent(in)                        :: AA                   !< beam cross section area
double precision, intent(in)                        :: I11                  !< beam cross section moments of interia about
double precision, intent(in)                        :: I22                  !< beam cross section moments of interia about
double precision, intent(in)                        :: Emat                 !< material youngs modulus
double precision, intent(in)                        :: btemp                !< beam temperature
double precision, intent(in)                        :: alpha                !< thermal expansion coefficient
!
! Output
!

!
! Input + Output
!
double precision, dimension(12), intent(inout)      :: disp_beam2           !< beam element displacement vector in global coordinates
double precision, dimension(12,12), intent(inout)   :: beam2_Kg             !< beam2-element geomtric stiffness matrix 
!
! Inner
!
double precision                                    :: facDP,le
! double precision                                    :: Ip
double precision, dimension(3)                      :: vv
double precision, dimension(12,12)                  :: TRMat
double precision                                    :: My1,My2,Mz1,Mz2,Fx,Fy,Fz
integer                                             :: jj,kk

integer                                             :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
facDP = 0.0
beam2_Kg = 0.0
!
! =================================================================================================
!
! Calculation

  ! get vector v (see MSC.NASTRAN 2005 Reference Guide p. 59)
  vv(1:3)=fesim%elemente%beam2s(ii)%xi(1:3)

  ! get transformation matrix from global to local coordinate system
  
  call beam2_rotation (beam2_node_coords,vv,'gl',TRMat,le)
  
  ! transform displacement vector to local coordinates
  
  disp_beam2=matmul(TRMat,disp_beam2)

  ! calculate element moments and forces
  
  !My1 = (Emat*I22/(le**2))*(6.D0*(disp_beam2(9)-disp_beam2(3))+4.D0*le*disp_beam2(5)+2.D0*le*disp_beam2(11))
  !My2 = (Emat*I22/(le**2))*(6.D0*(disp_beam2(9)-disp_beam2(3))+2.D0*le*disp_beam2(5)+4.D0*le*disp_beam2(11))
  
  !Mz1 = (Emat*I11/(le**2))*(-6.D0*(disp_beam2(8)-disp_beam2(2))+4.D0*le*disp_beam2(6)+2.D0*le*disp_beam2(12))
  !Mz2 = (Emat*I11/(le**2))*(-6.D0*(disp_beam2(8)-disp_beam2(2))+2.D0*le*disp_beam2(6)+4.D0*le*disp_beam2(12))
  
  My1 = 0.d0
  My2 = 0.d0
  
  Mz1 = 0.d0
  Mz2 = 0.d0
  
  Fx = (Emat*AA/le)*(disp_beam2(7)-disp_beam2(1)) - Emat*AA*btemp*alpha
  Fy = (-1.D0/le)*(Mz2+Mz1)
  Fz = ( 1.D0/le)*(My2+My1)
  
  ! polar moment of inertia

!   Ip = I11+I22
  
  ! fill lower diagonal of element stiffness matrix
  
  beam2_Kg(2,2)   =  6.D0*Fx/(5.D0*le)    ! 2nd row
  beam2_Kg(3,3)   =  beam2_Kg(2,2)        ! 3rd row
  beam2_Kg(4,2)   =  My1/le               ! 4th row
  beam2_Kg(4,3)   =  Mz1/le
  beam2_Kg(4,4)   =  0.d0 !Ip*Fx/(AA*le)
  beam2_Kg(5,3)   = -Fx/10.D0             ! 5th row
  beam2_Kg(5,4)   = -Fy*le/6.D0
  beam2_Kg(5,5)   =  2.D0*Fx*le/15.D0
  beam2_Kg(6,2)   = -beam2_Kg(5,3)        ! 6th row
  beam2_Kg(6,4)   = -Fz*le/6.D0
  beam2_Kg(6,6)   =  beam2_Kg(5,5)
  beam2_Kg(8,2)   = -beam2_Kg(2,2)        ! 8th row
  beam2_Kg(8,4)   = -beam2_Kg(4,2)
  beam2_Kg(8,6)   =  beam2_Kg(5,3)
  beam2_Kg(8,8)   =  beam2_Kg(2,2)
  beam2_Kg(9,3)   =  beam2_Kg(8,2)        ! 9th row
  beam2_Kg(9,4)   = -beam2_Kg(4,3)
  beam2_Kg(9,5)   =  beam2_Kg(6,2)
  beam2_Kg(9,9)   =  beam2_Kg(2,2)
  beam2_Kg(10,2)  =  My2/le               ! 10th row
  beam2_Kg(10,3)  =  Mz2/le
  beam2_Kg(10,4)  = -beam2_Kg(4,4)
  beam2_Kg(10,5)  = -beam2_Kg(5,4)
  beam2_Kg(10,6)  = -beam2_Kg(6,4)
  beam2_Kg(10,8)  = -beam2_Kg(10,2)
  beam2_Kg(10,9)  = -beam2_Kg(10,3)
  beam2_Kg(10,10) = -beam2_Kg(10,4)
  beam2_Kg(11,3)  =  beam2_Kg(8,6)        ! 11th row
  beam2_Kg(11,4)  = -beam2_Kg(5,4)
  beam2_Kg(11,5)  = -Fx*le/30.D0
  beam2_Kg(11,9)  =  beam2_Kg(6,2)
  beam2_Kg(11,10) =  beam2_Kg(5,4)
  beam2_Kg(11,11) =  beam2_Kg(5,5)
  beam2_Kg(12,2)  =  beam2_Kg(6,2)        ! 12th row
  beam2_Kg(12,4)  = -beam2_Kg(6,4)
  beam2_Kg(12,6)  =  beam2_Kg(11,5)
  beam2_Kg(12,8)  = -beam2_Kg(12,2)
  beam2_Kg(12,10) = -beam2_Kg(12,4)
  beam2_Kg(12,12) =  beam2_Kg(11,11)
  
  ! fill upper diagonal
  
  do jj=2,12
    do kk=(jj+1),12
      beam2_Kg(jj,kk) = beam2_Kg(kk,jj)
    end do
  end do
  
  ! transform element stiffness matrix to global coordinates
  ! Achtung die Matrix ist global nach lokal, deswegen das Transponiert gestauscht
  beam2_Kg = matmul(matmul(transpose(TRMat),beam2_Kg),TRMat)
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_geostiff'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine beam2_geostiff
