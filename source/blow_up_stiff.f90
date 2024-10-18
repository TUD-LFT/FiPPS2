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
!> Add zeros for rotational degrees in stiffness matrix of a solid element with n nodes
!
!> @details
!> Takes element stiffness matrix K_in with shape (3*n,3*n) of a solid element with
!> n nodes and blows it up to a (6*n,6*n) matrix with zeros at entries for rotational
!> degrees of freedom; by this, the blown up matrix can be handled with the degree of
!> freedom set to 6 for the following routines
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit 2015
!
!> $Id: blow_up_stiff.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine blow_up_stiff (n, K_in, K_out)

! =================================================================================================
! Use
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
!
integer, intent(in)                                :: n !< The element's number of nodes
double precision, dimension(3*n,3*n),intent(in)    :: K_in !< (3*n,3*n) element stiffness matrix of element with dof=3
!
! Output
!
double precision, dimension(6*n,6*n), intent(out)  :: K_out !< (6*n,6*n) element stiffness matrix of element, blown up to dof=6
!
! Internal
!
integer                                            :: ii, jj
integer                                            :: ii_s1, ii_e1, ii_s2, ii_e2
integer                                            :: jj_s1, jj_e1, jj_s2, jj_e2
!
! =================================================================================================
!
! Initialisation
!
K_out = 0.d0
!
! =================================================================================================
!
! Calculation
!
 do ii = 1,n
   do jj = 1,n
     ii_s1 = 6*(ii-1)+1
     ii_e1 = 6*(ii-1)+3
     ii_s2 = 3*(ii-1)+1
     ii_e2 = 3*(ii-1)+3
     jj_s1 = 6*(jj-1)+1
     jj_e1 = 6*(jj-1)+3
     jj_s2 = 3*(jj-1)+1
     jj_e2 = 3*(jj-1)+3
     K_out(ii_s1:ii_e1,jj_s1:jj_e1) = K_in(ii_s2:ii_e2,jj_s2:jj_e2)
   end do
 end do

end subroutine blow_up_stiff
