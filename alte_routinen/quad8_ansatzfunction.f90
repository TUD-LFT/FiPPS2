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
!> Calculates all necessary values for numerical integration of a quad8 finite
!> shell element
!
!> @details
!
!> @author   Martin Rädel, TU Dresden, Diplomarbeit, 08.12.2009
!               
!> @author   Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 17.10.2012
!
! =================================================================================================
subroutine quad8_ansatzfunction(xi, eta, node_coords, withNi, Ni, dNidx, dNidy, detJac)
! =================================================================================================
!
! Use
!

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
    subroutine quad8_ansatzfunction_xieta(xi, eta, Ni, dNidxi, dNideta, d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2)
      double precision, intent(in)                              :: xi, eta
      
      double precision, dimension(8), optional, intent(out)     :: Ni
      double precision, dimension(8), optional, intent(out)     :: dNidxi, dNideta
      double precision, dimension(8), optional, intent(out)     :: d2Nidxi2, d2Nidxideta, d2Nidetadxi, d2Nideta2
    end subroutine quad8_ansatzfunction_xieta
  end interface
!
! =================================================================================================
!
! Data types
!
  double precision, intent(in)                                  :: xi           !< xi coordinate inside the element (-1 <= xi <= 1)
  double precision, intent(in)                                  :: eta          !< eta coordinate inside the element (-1 <= eta <= 1)
  double precision, dimension(8,4), intent(in)                  :: node_coords  !< (8,4)-array containing the coordinates of each node (
  logical, intent(in)                                           :: withNi       !< flag if Ni should be calculated or not
  double precision, dimension(8), intent(out)                   :: Ni           !< (8)-array containing the values of the shape function for each node
  double precision, dimension(8), intent(out)                   :: dNidx        !< value of the first derivation of the shape function which respect to x for each node
  double precision, dimension(8), intent(out)                   :: dNidy        !< value of the first derivation of the shape function which respect to y for each node
  double precision, intent(out)                                 :: detJac       !< Determinant of the jacobi matrix
  
  integer                                                       :: jj

  double precision                                              :: facDP, dxdxi, dxdeta, dydxi, dydeta
  double precision, dimension(8)                                :: dNidxi,dNideta
  double precision                                              :: dxidx, dxidy, detadx,detady
  
  integer                                                       :: err_code=0
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
  ! calculate shape-function values with respect to natural coordinates
  
  if (withNi .eq. .true.) then 
    
    call quad8_ansatzfunction_xieta(xi, eta, Ni=Ni, dNidxi=dNidxi, dNideta=dNideta)
    
  else
  
    call quad8_ansatzfunction_xieta(xi, eta, dNidxi=dNidxi, dNideta=dNideta)
  
  end if
  
  ! Calculate Jacobian
  
  dxdxi  = 0.D0
  dxdeta = 0.D0
  dydxi  = 0.D0
  dydeta = 0.D0
  
  do jj=1,8
     
    dxdxi  = dxdxi +dNidxi(jj) *node_coords(jj,2)
    dxdeta = dxdeta+dNideta(jj)*node_coords(jj,2)
    
    dydxi  = dydxi +dNidxi(jj) *node_coords(jj,3)
    dydeta = dydeta+dNideta(jj)*node_coords(jj,3)
     
  end do
  
  detJac = dxdxi*dydeta-dxdeta*dydxi
  
  facDP = 1.D0/detJac
 
  dxidx  =  facDP*dydeta
  dxidy  = -facDP*dxdeta
  detadx = -facDP*dydxi
  detady =  facDP*dxdxi
  
  ! Calculate shape function derivatives with respect to local cartesian coordinate system
  
  do jj=1,8
     
    dNidx(jj) = dNidxi(jj)*dxidx+dNideta(jj)*detadx   ! == delta1
    dNidy(jj) = dNidxi(jj)*dxidy+dNideta(jj)*detady   ! == delta2
    
  end do

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_ansatzfunction'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine quad8_ansatzfunction
