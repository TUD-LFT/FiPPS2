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
!
!> @details
!
!> @author 
!
!> $Id: vtkouteig.f90 362 2018-08-07 09:18:14Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 362 $
!> $Date: 2018-08-07 11:18:14 +0200 (Di, 07. Aug 2018) $
!
! =================================================================================================
SUBROUTINE vtkoutelemcoord(fesim)
!
! Use
!
  use fesimulation_typen
  USE vtk_variablen
!
! =================================================================================================
!
  IMPLICIT NONE
!
! =================================================================================================
!
! Interface
!

!
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
!
! Input
!
  type(fe_simulation)                                   :: fesim
!
! Internal
!
  integer                                               :: ii, jj, cid, offset
  
  integer                                               :: err_code=0
  
  double precision, dimension(:,:), allocatable         :: xVec
  double precision, dimension(:,:), allocatable         :: yVec
  double precision, dimension(:,:), allocatable         :: zVec
  double precision, dimension(3)                        :: xElem
  double precision, dimension(3)                        :: yElem
  double precision, dimension(3)                        :: vv
  double precision, dimension(2,3)                      :: beam2_node_coords
  double precision, dimension(20,3)                     :: lsolid20_node_coords
  double precision, dimension(20)                       :: dNidx,dNidy,dNidz
  double precision, dimension(3,3)                      :: J, E
  double precision                                      :: detJac,absVec,le
!
! =================================================================================================
!
! Initialisation
!
  offset = 0
  
  allocate(xVec(fesim%num_elements,3))
  allocate(yVec(fesim%num_elements,3))
  allocate(zVec(fesim%num_elements,3))
!
! =================================================================================================
!
! Calculation
!  
  xVec = 0.d0
  yVec = 0.d0
  zVec = 0.d0
  
  IF (ALLOCATED(fesim%elemente%beam2s) .EQ. .TRUE.) THEN
    DO ii = 1, size(fesim%elemente%beam2s,1)

      DO jj = 1,2
        beam2_node_coords(jj,:) = fesim%knoten%nodes(fesim%elemente%beam2s(ii)%nids(jj))%coords(1:3)
      END DO

      vv(1:3)=fesim%elemente%beam2s(ii)%xi(1:3)

      ! Transformation matrix from element to global coordinate system
      CALL beam2_cosinematrix(beam2_node_coords,vv,E,le)

      xVec(offset + ii,1:3) = E(:,1)
      yVec(offset + ii,1:3) = E(:,2)
      zVec(offset + ii,1:3) = E(:,3)

    END DO
    offset = offset + size(fesim%elemente%beam2s,1)
  END IF
  
  IF (ALLOCATED(fesim%elemente%quad8s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%quad8s,1)
      
      xElem(1) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(2))%coords(1) - fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(1))%coords(1)
      xElem(2) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(2))%coords(2) - fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(1))%coords(2)
      xElem(3) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(2))%coords(3) - fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(1))%coords(3)
      
      absVec = SQRT(xElem(1)**2 + xElem(2)**2 + xElem(3)**2)
      
      xElem = xElem/absVec
      
      yElem(1) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(4))%coords(1) - fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(1))%coords(1)
      yElem(2) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(4))%coords(2) - fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(1))%coords(2)
      yElem(3) = fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(4))%coords(3) - fesim%knoten%nodes(fesim%elemente%quad8s(ii)%nids(1))%coords(3)
      
      absVec = SQRT(yElem(1)**2 + yElem(2)**2 + yElem(3)**2)
      
      yElem = yElem/absVec
    
      xVec(offset + ii,1) = xElem(1)
      xVec(offset + ii,2) = xElem(2)
      xVec(offset + ii,3) = xElem(3)
    
      zVec(offset + ii,1) = xElem(2)*yElem(3) - xElem(3)*yElem(2)
      zVec(offset + ii,2) = xElem(3)*yElem(1) - xElem(1)*yElem(3)
      zVec(offset + ii,3) = xElem(1)*yElem(2) - xElem(2)*yElem(1)
      
      yVec(offset + ii,1) = zVec(offset + ii,2)*xElem(3) - zVec(offset + ii,3)*xElem(2)
      yVec(offset + ii,2) = zVec(offset + ii,3)*xElem(1) - zVec(offset + ii,1)*xElem(3)
      yVec(offset + ii,3) = zVec(offset + ii,1)*xElem(2) - zVec(offset + ii,2)*xElem(1)
      
    END DO
    offset = offset + size(fesim%elemente%quad8s,1)
  END IF
    
  IF (ALLOCATED(fesim%elemente%lsolid20s) .EQ. .TRUE.) THEN 
    DO ii = 1, size(fesim%elemente%lsolid20s,1)

      DO jj = 1,20
        lsolid20_node_coords(jj,:) = fesim%knoten%nodes(fesim%elemente%lsolid20s(ii)%nids(jj))%coords(1:3)
      END DO

      ! Get derivates of shape functions in element center with respect to global coordinates
      CALL hexa20_shapefunctions_x_y_z(0.d0,0.d0,0.d0,lsolid20_node_coords,dNidx,dNidy,dNidz,detJac,J)
 
      ! get ID of coordinate system for orientation of element coosy
      cid  = fesim%eigenschaften%plsolids(fesim%elemente%lsolid20s(ii)%int_pid)%cid

      ! Transformation matrix from global to element coordinate system
      CALL cosinematrix_esys(fesim,J,cid,E)

      xVec(offset + ii,1:3) = E(1,:)
      yVec(offset + ii,1:3) = E(2,:)
      zVec(offset + ii,1:3) = E(3,:)

    END DO

    offset = offset + size(fesim%elemente%lsolid20s,1)
  END IF
  
  CALL vtkoutcelldata_vector(xVec, fesim%num_elements, 'ElemCoord_X')
  CALL vtkoutcelldata_vector(yVec, fesim%num_elements, 'ElemCoord_Y')
  CALL vtkoutcelldata_vector(zVec, fesim%num_elements, 'ElemCoord_Z')
  
  deallocate(xVec, yVec, zVec)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'vtkoutelemcoord'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE vtkoutelemcoord
