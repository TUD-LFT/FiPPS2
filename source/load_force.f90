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
!> subroutine for computing load-vectors entrys resulting from single forces
!
!> @details
!> subroutine calculates the force contribution to the loadvector
!
!> @author Martin Rädel, TU Dresden, Diplomarbeit, 08.12.2009
!
!> @author Florian Dexl, TU Dresden, Diplomarbeit, 09.10.2015
!>                                  überflüssige Schleife über nodes(:) entfernt,
!>                                  zur Einsparung von Rechenzeit
!
!> $Id: load_force.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine load_force (fesim,hh,lvec,pos)


use konstanten
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
type(fe_simulation)                         :: fesim
integer                                     :: pos_begin, pos_end, pos

double precision, dimension(:,:)            :: lvec(1:fesim%num_dof,1:fesim%num_loadcases)
double precision, dimension(:)              :: load_i(3)

integer                                     :: hh,ii
integer                                     :: err_code=0
!
! =================================================================================================
!
! Initialisation
!

!
! =================================================================================================
!
! Calculation

do ii=1,size(fesim%lasten%forces,1)                                             ! Loop over forces
   if (fesim%lasten%loads(hh)%lidi == fesim%lasten%forces(ii)%lid) then
   
      load_i(1:3)=fesim%lasten%forces(ii)%fac*fesim%lasten%forces(ii)%ni(1:3)
            
      ! Scale load with factor sfaci for single load in load card
      
      load_i = fesim%lasten%loads(hh)%sfaci*load_i
      
      ! Scale load with factor sfac for whole loadcase in load card
      
      load_i = fesim%lasten%loads(hh)%sfac*load_i

      pos_begin=(fesim%lasten%forces(ii)%nid-1)*ndof+1
      pos_end  =pos_begin+2         

      lvec(pos_begin:pos_end,pos)=lvec(pos_begin:pos_end,pos)+load_i

   end if
end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'load_force'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine load_force
