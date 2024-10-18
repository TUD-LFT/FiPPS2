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
!> Updates the geometry with the displacements of the previous calculation
!
!> @details
!> The translational displacements u_x, u_y, u_z of the previous calculation are
!> added to the nodal coordinates; performs the update of the geometry using the
!> displacements of the pervious calculation; this is just done, if the logical
!> "upgeom" for the current step of a multistep-calculation is set true
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: update_geometry.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine update_geometry(fesim, Utot)
!
! use
!
use konstanten
use fesimulation_typen
!
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
type(fe_simulation)                                     :: fesim
double precision, dimension(fesim%num_dof), intent(in)  :: Utot !< Array containing the displacements of the previous calculation at all nodes
!
! Internal
!
integer                                                 :: err_code=0, ii, n_1, n_2
!
! =================================================================================================
!
! Initialisation
!
! =================================================================================================
!
! Calculation
!
 do ii = 1, fesim%num_nodes

   n_1 = (ii-1)*ndof+1
   n_2 = n_1+2
 
   fesim%knoten%nodes(ii)%coords(1:3) = fesim%knoten%nodes(ii)%coords(1:3) + Utot(n_1:n_2)
 
 end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'update_geometry'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine update_geometry
