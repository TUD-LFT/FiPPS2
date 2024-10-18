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
!> Shape functions of the beam2 element for diplacements.
!
!> @details 
!> Calculates the shape functions of the beam2 element.
!> For bending displacements in elemental y- and z-direction, a cubic shape function is used.
!> For elongation in elemental x-direction, a linear shape function is used, as well as for the
!> torsional rotation. The rotations about elemental y- and z-directions are calculated as the
!> derivates of the bending displacements. Therefore, also the derivatives of the shape functions
!> with respect to xi can be requested.
!
!> @author Florian Dexl, TU Dresden, wissenschaftlicher Mitarbeiter
!
!> @date
!> 17.10.2019
!
!> $Id: beam2_rotation.f90 362 2018-08-07 09:18:14Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 362 $
!> $Date: 2018-08-07 11:18:14 +0200 (Di, 07. Aug 2018) $
! 
! =================================================================================================
subroutine beam2_shape_functions (xi,le,Ncubic1,Ncubic2,Nlinear,dNcubic1dxi,dNcubic2dxi,dNlineardxi)
!
use konstanten
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
double precision, intent(in)                          :: xi                  !< normalized running coordinate of beam element
double precision, intent(in)                          :: le                  !< length of beam element
!
! Output
!
double precision, dimension(2), optional, intent(out) :: Ncubic1             !< values of the first cubic shape function
double precision, dimension(2), optional, intent(out) :: Ncubic2             !< values of the second cubic shape function
double precision, dimension(2), optional, intent(out) :: Nlinear             !< values of the linear shape function
double precision, dimension(2), optional, intent(out) :: dNcubic1dxi         !< values of the first derivation of the first cubic shape function which respect to xi
double precision, dimension(2), optional, intent(out) :: dNcubic2dxi         !< values of the first derivation of the second cubic shape function which respect to xi
double precision, dimension(2), optional, intent(out) :: dNlineardxi         !< values of the first derivation of the linear shape function which respect to xi
!
! Internal
!
integer                                          :: err_code=0
!
! =================================================================================================
!
! Initialisation
!
! =================================================================================================
!
! Calculation
!
  ! calculation of the first cubic shape function for bending displacement
  if (present(Ncubic1) .eq. .true.) then
    Ncubic1(1)  =   2.d0*xi**3 - 3.d0*xi**2 +  1
    Ncubic1(2)  =  -2.d0*xi**3 + 3.d0*xi**2
  end if

  ! calculation of the second cubic shape function for bending displacement
  if (present(Ncubic2) .eq. .true.) then
    Ncubic2(1)  = (      xi**3 - 2.d0*xi**2 + xi)*le
    Ncubic2(2)  = (      xi**3 -      xi**2     )*le
  end if

  ! calculation of linear shape functions for elongation and torsional rotation
  if (present(Nlinear) .eq. .true.) then
    Nlinear(1) = 1.d0 - xi
    Nlinear(2) = xi
  end if
  
  if (present(dNcubic1dxi) .eq. .true.) then
    dNcubic1dxi(1) =   6.d0*xi**2 - 6.d0*xi
    dNcubic1dxi(2) =  -6.d0*xi**2 + 6.d0*xi
  end if
  
  if (present(dNcubic2dxi) .eq. .true.) then
    dNcubic2dxi(1) = ( 3.d0*xi**2 - 4.d0*xi + 1)*le
    dNcubic2dxi(2) = ( 3.d0*xi**2 - 2.d0*xi    )*le
  end if

  if (present(dNlineardxi) .eq. .true.) then
    dNlineardxi(1) = -1.d0
    dNlineardxi(2) =  1.d0
  end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_shape_functions'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine beam2_shape_functions
