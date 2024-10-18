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
!> Calculates coupling matrix for interpolation of the degrees of freedom at the beam2 element
!
!> @details
!> Calculates coupling matrix for interpolation of the degrees of freedom at the beam2 element.
!> For the given value xi (normalized running coordinate along beam2 element), a coupling matrix
!> is calculated that allows the interpolation of the degrees of freedom at the point xi out
!> of the degrees of freedom of the two nodes of the beam2 element.
!> The coupling matrix is given with respect to the elemental coordinate system.
!
!> @author
!> Florian Dexl, TU Dresden, wissenschaftlicher Mitarbeiter
!
!> @date
!> 21.10.2019
!
!> $Id: impc_contact_node_beam2.f90 383 2018-11-02 08:07:44Z DOM+ahauffe $
!> $Author: DOM+ahauffe $
!> $Revision: 383 $
!> $Date: 2018-11-02 09:07:44 +0100 (Fr, 02. Nov 2018) $
!
! =================================================================================================
subroutine beam2_couplingmatrix (xi, le, couplingMatrix, Nlinear)
!
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Interface
!
  interface 
    subroutine beam2_shape_functions(xi, le, Ncubic1, Ncubic2, Nlinear, dNcubic1dxi, dNcubic2dxi, dNlineardxi)
      double precision, intent(in)                              :: xi, le
      
      double precision, dimension(2), optional, intent(out)     :: Ncubic1
      double precision, dimension(2), optional, intent(out)     :: Ncubic2
      double precision, dimension(2), optional, intent(out)     :: Nlinear
      double precision, dimension(2), optional, intent(out)     :: dNcubic1dxi
      double precision, dimension(2), optional, intent(out)     :: dNcubic2dxi
      double precision, dimension(2), optional, intent(out)     :: dNlineardxi
    end subroutine beam2_shape_functions
  end interface
!
! =================================================================================================
!
! Input
!
double precision, intent(in)                         :: xi, le
!
! Output
!
double precision, dimension(6,12), intent(out)       :: couplingMatrix
double precision, dimension(2), intent(out), optional:: Nlinear
!
! inner
!
integer                                              :: ii
double precision, dimension(2)                       :: Ncubic1, Ncubic2
double precision, dimension(2)                       :: dNcubic1dxi, dNcubic2dxi
double precision, dimension(2)                       :: Nlinear_
integer                                              :: err_code=0

  call beam2_shape_functions(xi, le, Ncubic1=Ncubic1, Ncubic2=Ncubic2, Nlinear=Nlinear_, dNcubic1dxi=dNcubic1dxi, dNcubic2dxi=dNcubic2dxi)

  ! Zusammenstellen der Kopplungsmatrix
  couplingMatrix(:,:) = 0.d0
  do ii = 1,2
    couplingMatrix(1,(ii-1)*6+1) = Nlinear_(ii)
    couplingMatrix(2,(ii-1)*6+2) = Ncubic1(ii)
    couplingMatrix(2,(ii-1)*6+6) = Ncubic2(ii)
    couplingMatrix(3,(ii-1)*6+3) = Ncubic1(ii)
    couplingMatrix(3,(ii-1)*6+5) = -1.d0*Ncubic2(ii)
    couplingMatrix(4,(ii-1)*6+4) = Nlinear_(ii)
    couplingMatrix(5,(ii-1)*6+3) = -1.d0*dNcubic1dxi(ii)/le
    couplingMatrix(5,(ii-1)*6+5) = dNcubic2dxi(ii)/le
    couplingMatrix(6,(ii-1)*6+2) = dNcubic1dxi(ii)/le
    couplingMatrix(6,(ii-1)*6+6) = dNcubic2dxi(ii)/le
  end do

  if (present(Nlinear) .eq. .true.) then
    Nlinear = Nlinear_
  end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'beam2_couplingmatrix'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine beam2_couplingmatrix
