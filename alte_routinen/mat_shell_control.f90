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
subroutine mat_shell_control (jj,propid,flag,mattype,propm1,Amat,Bmat,Dmat,Hmat,SOmega)
! =================================================================================================
!
!	Header:		control subroutine for getting material properties
!
!	Content:	
!
!	Input:		
!			flag:	flag = true:	get Amat,Bmat,Dmat,Hmat,SOmega for element
!						stiffness matrix
! 				flag = false:	get Amat and Bmat only for element geometric
!						stiffness matrix
!
!	Output:		
!
!	Internal:	
!
!	Calls:		
!
!	Called by:	quad8_stiff_control,
!			quad8_geostiff_control
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
!	Revision:	Martin Rädel			02.07.2010
!			TU Dresden, WiMi
!			-> flag eingeführt
!
! =================================================================================================
!
! Use
!
use Globale_Variablen
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
integer, intent(in)				:: jj,propid
logical, intent(in)				:: flag
!
! Output
!
double precision, dimension(3,3), intent(out)	:: Amat,Bmat,Dmat
double precision, dimension(2,2), intent(out)	:: Hmat
double precision, intent(out)			:: SOmega
!
! Input + Output
!
 character(3), intent(inout)			:: mattype
integer, intent(inout)				:: propm1
!
! inner
!
integer						:: kk
integer						:: err_code=0
double precision                                :: thick
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
   if (propid /= propm1) then
         
      ! isotropic material
      
      if (is_pshell == .true.) then
         do kk=1,size(pshells,1)
	    
	    if (propid == pshells(kk)%pid) then
               
               call mat_shell_iso (flag,kk,Amat,Bmat,Dmat,Hmat,SOmega)
	       
	       mattype='iso'
               
            end if
         end do
      end if
      
      ! layerwise othotropic material
      
      if (is_pcomp == .true.) then
         do kk=1,size(pcomps,1)

	    if (propid == pcomps(kk)%pid) then
               
               call mat_shell_ortho (.true.,kk,Amat,Bmat,Dmat,Hmat,SOmega,thick)
               pcomps(kk)%thick = thick
	       
	       mattype='ort'
               
            end if
         end do
      end if
      
      propm1 = propid
      
   end if
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)  			   'mat_shell_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine mat_shell_control
