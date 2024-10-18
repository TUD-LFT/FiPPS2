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
!> Calculates the strains for an 8-node shell element quad8 at all integration points
!
!> @details
!> Calculation of the strain for an 8-node shell element quad8 at all integration 
!> points at the position zeta in thickness direction; strains are in the element
!> coordinate system
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.05.2017
!
!> $Id: quad8_ip_strain.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine quad8_ip_strain (node_coords, disp_quad8, zeta, nip, tth, area, C31, C32, C33, C34, ntemps, ath, strain, is_temperature, is_quad8temp)
! =================================================================================================
! use
!
  use integration_schemes
  
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
! Input
!
  double precision, intent(in)                      :: node_coords(8,3)   !< (8,3)-Array containing global coordinates of element nodes
  double precision, intent(in)                      :: disp_quad8(48)     !< (48)-Array containing als resultant displacemants in order (u1, v1, w1, rotx1, roty1, rotz1, u2, ...)
  double precision, intent(in)                      :: zeta               !< relative position in thickness direction
  integer, intent(in)                               :: nip                !< Number of integration points in plane per direction
  double precision, intent(in)                      :: tth                !< Total thickness of laminat (sum of all values in lth)
  double precision, intent(in)                      :: area               !< element reference surface area
  double precision, intent(in)                      :: C31                !< 31-term of the 3D material stiffness matrix C to calculate the strain in normal direction
  double precision, intent(in)                      :: C32                !< 32-term of the 3D material stiffness matrix C to calculate the strain in normal direction
  double precision, intent(in)                      :: C33                !< 33-term of the 3D material stiffness matrix C to calculate the strain in normal direction
  double precision, intent(in)                      :: C34                !< 34-term of the 3D material stiffness matrix C to calculate the strain in normal direction
  double precision, intent(in)                      :: ntemps(8)          !< (8)-Array containing the temperature at each node
  double precision, intent(in)                      :: ath(6)             !< (6)-Array containing the coefficients of thermal expansion
  logical, intent(in)                               :: is_temperature
  logical, intent(in)                               :: is_quad8temp
!
! Output
!
  double precision, dimension(nip*nip,6), intent(out)   :: strain         !< (nip*nip,6)-Array containing the strain tensor at each inplane integration point at relative thickness position zeta. The shear strain are gamma = 2 * epsilon

  double precision, dimension(8)                    :: Ni
  double precision, dimension(3,8)                  :: normalVec                ! Normalenvektoren an den acht Punkten
  double precision, dimension(8,3,3)                :: nodeTransMat             ! lokales Koordinatensystem an jedem Knoten
  
  double precision, dimension(:), allocatable       :: xi_vec, eta_vec, w_xi, w_eta
  
  double precision, dimension(3,3,8)                :: nodeTransMat2
  
  double precision, dimension(48,6)                 :: bMat
  
  double precision, dimension(3,3)                  :: Jac
  
  double precision, dimension(3)                    :: xElem
  
  double precision                                  :: k, detJac
  
  integer                                           :: ii, jj
  integer                                           :: numIntPointsIP
  double precision, dimension(3,3)                  :: transMat
  double precision, dimension(:)                    :: reorderedDispl(48)
  double precision, dimension(6)                    :: eps_th
  
  integer                                           :: err_code=0
  
  do jj = 1,6
    do ii = 1,8
      reorderedDispl((jj-1)*8+ii) = disp_quad8((ii-1)*6+jj)
    end do
  end do
  
  ! Berechnung der Normalen auf der Elementmittelfläche und der dazugehörigen Koordinatensysteme mit Transformationsmatrizen
  call quad8_normals(node_coords,normalVec,nodeTransMat,nodeTransMat2)
    
  xElem(1) = node_coords(2,1) - node_coords(1,1)
  xElem(2) = node_coords(2,2) - node_coords(1,2)
  xElem(3) = node_coords(2,3) - node_coords(1,3)
  
  numIntPointsIP = nip*nip
    
  allocate(xi_vec(numIntPointsIP), eta_vec(numIntPointsIP))
  allocate(  w_xi(numIntPointsIP),   w_eta(numIntPointsIP))
  
  call integration_points_2d(nip, "Gauss",  xi_vec, w_xi, &
                             nip, "Gauss", eta_vec, w_eta)
  
  k = MAX(1.2d0, 1.d0+0.2d0*area/(25.d0*tth**2))
  
  do ii = 1,numIntPointsIP
  
    ! Berechnung des Koordinatensystems im Integrationspunkt (nur in der Mittelebene, gilt für die gesamte Dickenrichtung)
    call quad8_ip_coordsys(node_coords, xi_vec(ii), eta_vec(ii), xElem, transMat)
    
    call quad8_ansatzfunction_xieta(xi_vec(ii), eta_vec(ii), Ni=Ni)
            
    ! Beschleunigte Berechung der B-Matrix durch Verzicht auf viele der Nulloperationen
    call quad8_bmat(xi_vec(ii), eta_vec(ii), zeta, node_coords, tth, normalVec, nodeTransMat, bmat, detJac, Jac, transMat)

    strain(ii,1:6) = MatMul(transpose(Bmat),reorderedDispl)
    
    ! epsilon ist über die Bedingung 0 = simga33 = C13 * epsilon11 + C23 * epsilon22 + C33 * epsilon33 aus 
    ! dem Materialgesetz und dem ebenen Spannungszustand hergeleitet
    strain(ii,3)   = - C31/C33 * strain(ii,1) - C32/C33 * strain(ii,2) - C34/C33 * strain(ii,4)
    
    strain(ii,5)   = strain(ii,5)/k
    strain(ii,6)   = strain(ii,6)/k
    
    if ((is_temperature == .true.) .or. (is_quad8temp == .true.)) then
      eps_th(1:6)  = SUM(ntemps(1:8)*Ni(1:8)) * ath(1:6)
      strain(ii,3) = - C31/C33 * (strain(ii,1) - eps_th(1)) - C32/C33 * (strain(ii,2) - eps_th(2)) - C34/C33 * (strain(ii,4) - eps_th(4)) 
      strain(ii,1:6) = strain(ii,1:6) - eps_th(1:6)
    end if
  
  end do

  deallocate(xi_vec, eta_vec)
  deallocate(  w_xi,   w_eta)

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'quad8_ip_strain_stress'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine quad8_ip_strain
