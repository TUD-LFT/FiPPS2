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
!> $Id: nodalCoordSysTransMat.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine nodalCoordSysTransMat(fesim,num_node,node_ids,transMat)

  use konstanten
  use fesimulation_typen

  implicit none
  
  type(fe_simulation)                                             :: fesim
  integer                                                         :: num_node
  integer                                                         :: node
  integer, dimension(1:num_node)                                  :: node_ids
  
  double precision, dimension(num_node*ndof,num_node*ndof)        :: transMat
  
  transMat = 0.d0
  do node = 1, num_node
 
    if (fesim%knoten%nodes(node_ids(node))%cid .eq. 0) then
      transMat((node-1)*ndof+1,(node-1)*ndof+1) = 1.d0
      transMat((node-1)*ndof+2,(node-1)*ndof+2) = 1.d0
      transMat((node-1)*ndof+3,(node-1)*ndof+3) = 1.d0
      transMat((node-1)*ndof+4,(node-1)*ndof+4) = 1.d0
      transMat((node-1)*ndof+5,(node-1)*ndof+5) = 1.d0
      transMat((node-1)*ndof+6,(node-1)*ndof+6) = 1.d0
    else      
      transMat((node-1)*ndof+1:(node-1)*ndof+3,(node-1)*ndof+1:(node-1)*ndof+3) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(node_ids(node))%cid)%transMat(1:3,1:3)
      transMat((node-1)*ndof+4:(node-1)*ndof+6,(node-1)*ndof+4:(node-1)*ndof+6) = fesim%koordinatensysteme%coords(fesim%knoten%nodes(node_ids(node))%cid)%transMat(1:3,1:3)
    end if
  
  end do
  
end subroutine nodalCoordSysTransMat
