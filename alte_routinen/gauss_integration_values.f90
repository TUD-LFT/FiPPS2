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
subroutine gauss_integration_values (numPoints, xi_vec, eta_vec, w_xi, w_eta)
! =================================================================================================
!
!	Header:		Liefert alle notwendigen Daten für eine Gausspunktintegration
!
!	Content:	Gibt die Koordinaten der Stützenstellen und die dazugehörigen Gewichte
!                       für eine 1,2 oder 3- Gausspunktintragration zurück
!
!	Input:		
!
!	Output:		
!
!	Internal:	
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!			Andreas Hauffe			16.06.2010
! 			TU Dresden, WiMi
!
!	Revision:	
!
! =================================================================================================
!
! use
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
  integer, intent(in)                                           :: numPoints
  integer					                :: err_code = 0
  
  double precision				                :: facDP
  double precision, dimension(numPoints*numPoints), intent(out) :: xi_vec, eta_vec
  double precision, dimension(numPoints*numPoints), intent(out) :: w_xi, w_eta
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
  if (numPoints == 3) then
     
    ! Build vectors with Gauß-point coordinates in xi and eta
    
    facDP = sqrt(0.6D0)
    
    xi_vec(1) = -facDP
    xi_vec(2) =  0.D0
    xi_vec(3) =  facDP
    xi_vec(4) = -facDP
    xi_vec(5) =  0.D0
    xi_vec(6) =  facDP
    xi_vec(7) = -facDP
    xi_vec(8) =  0.D0
    xi_vec(9) =  facDP
    
    eta_vec(1) = -facDP
    eta_vec(2) = -facDP
    eta_vec(3) = -facDP
    eta_vec(4) =  0.D0
    eta_vec(5) =  0.D0
    eta_vec(6) =  0.D0
    eta_vec(7) =  facDP
    eta_vec(8) =  facDP
    eta_vec(9) =  facDP
    
    ! Build vectors with Gauss-integration weights
    
    w_xi(1) = 5.D0/9.D0
    w_xi(2) = 8.D0/9.D0
    w_xi(3) = w_xi(1)
    w_xi(4) = w_xi(1)
    w_xi(5) = w_xi(2)
    w_xi(6) = w_xi(1)
    w_xi(7) = w_xi(1)
    w_xi(8) = w_xi(2)
    w_xi(9) = w_xi(1)
    
    w_eta(1) = 5.D0/9.D0
    w_eta(2) = w_eta(1)
    w_eta(3) = w_eta(1)
    w_eta(4) = 8.D0/9.D0
    w_eta(5) = w_eta(4)
    w_eta(6) = w_eta(4)
    w_eta(7) = w_eta(1)
    w_eta(8) = w_eta(1)
    w_eta(9) = w_eta(1)
  
  else if (numPoints == 2) then
     
    facDP = 1.D0/sqrt(3.D0)
    
    xi_vec(1) = -facDP
    xi_vec(2) =  facDP
    xi_vec(3) =  facDP
    xi_vec(4) = -facDP
    
    eta_vec(1) = -facDP
    eta_vec(2) = -facDP
    eta_vec(3) =  facDP
    eta_vec(4) =  facDP
    
    w_xi(1)  = 1.D0
    w_xi(2)  = w_xi(1)
    w_xi(3)  = w_xi(1)
    w_xi(4)  = w_xi(1)
    
    w_eta(1) = 1.D0
    w_eta(2) = w_eta(1)
    w_eta(3) = w_eta(1)
    w_eta(4) = w_eta(1)
     
  else if (numPoints == 1) then
     
    facDP      = 0.D0
    
    xi_vec(1)  = facDP
    
    eta_vec(1) = facDP
    
    w_xi(1)    = 2.D0
    
    w_eta(1)   = 2.D0
     
  else
     
    write(*,*) 'wrong input on parameter numPoints'
    write(*,*) 'values 1, 2 or 3 are allowed to be chosen'
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
    write(*,*)  		    'gauss_integration_values'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine gauss_integration_values
