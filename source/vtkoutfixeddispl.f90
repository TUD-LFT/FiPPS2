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
SUBROUTINE vtkoutfixeddispl(fesim)
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
  integer                                               :: ii
  
  integer                                               :: err_code=0
  
  double precision, dimension(:,:), allocatable         :: xDis,yDis,zDis
  double precision, dimension(:,:), allocatable         :: xRot,yRot,zRot
  logical                                               :: fixedNode
!
! =================================================================================================
!
! Initialisation
!
  allocate(xDis(fesim%num_nodes,3))
  allocate(yDis(fesim%num_nodes,3))
  allocate(zDis(fesim%num_nodes,3))
  allocate(xRot(fesim%num_nodes,3))
  allocate(yRot(fesim%num_nodes,3))
  allocate(zRot(fesim%num_nodes,3))
!
! =================================================================================================
!
! Calculation
!  

  xdis = 0.d0
  ydis = 0.d0
  zdis = 0.d0
  xrot = 0.d0
  yrot = 0.d0
  zrot = 0.d0
  
  DO ii = 1,fesim%num_nodes
    fixedNode = .FALSE.
    IF (fesim%internals%num_dof_vec(6*(ii-1)+1) .EQ. 0) THEN
      xdis(ii,1) = 1.d0
    fixedNode = .TRUE.
    END IF
    IF (fesim%internals%num_dof_vec(6*(ii-1)+2) .EQ. 0) THEN
      ydis(ii,2) = 1.d0
    fixedNode = .TRUE.
    END IF
    IF (fesim%internals%num_dof_vec(6*(ii-1)+3) .EQ. 0) THEN
      zdis(ii,3) = 1.d0
    fixedNode = .TRUE.
    END IF
    IF (fesim%internals%num_dof_vec(6*(ii-1)+4) .EQ. 0) THEN
      xrot(ii,1) = 1.d0
    fixedNode = .TRUE.
    END IF
    IF (fesim%internals%num_dof_vec(6*(ii-1)+5) .EQ. 0) THEN
      yrot(ii,2) = 1.d0
    fixedNode = .TRUE.
    END IF
    IF (fesim%internals%num_dof_vec(6*(ii-1)+6) .EQ. 0) THEN
      zrot(ii,3) = 1.d0
    fixedNode = .TRUE.
    END IF
    
    if (fixedNode .EQ. .TRUE. .AND. fesim%knoten%nodes(ii)%cid .ne. 0) then
      xdis(ii,1:3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),xdis(ii,1:3))
      ydis(ii,1:3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),ydis(ii,1:3))
      zdis(ii,1:3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),zdis(ii,1:3))
      xrot(ii,1:3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),xrot(ii,1:3))
      yrot(ii,1:3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),yrot(ii,1:3))
      zrot(ii,1:3) = matmul(fesim%koordinatensysteme%coords(fesim%knoten%nodes(ii)%cid)%transMat(1:3,1:3),zrot(ii,1:3))
    end if
    
  END DO

  CALL vtkoutpointvektor(xdis, fesim%num_nodes, 'FixedDis_X')
  CALL vtkoutpointvektor(ydis, fesim%num_nodes, 'FixedDis_Y')
  CALL vtkoutpointvektor(zdis, fesim%num_nodes, 'FixedDis_Z')
  CALL vtkoutpointvektor(xrot, fesim%num_nodes, 'FixedRot_X')
  CALL vtkoutpointvektor(yrot, fesim%num_nodes, 'FixedRot_Y')
  CALL vtkoutpointvektor(zrot, fesim%num_nodes, 'FixedRot_Z')
  
  deallocate(xDis,yDis,zDis,xRot,yRot,zRot)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'vtkoutfixeddispl'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

END SUBROUTINE vtkoutfixeddispl
