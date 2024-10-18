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
subroutine transform_nodal_coords(tMat, num_dof_vec, dim_dof)
! =================================================================================================
!
!	Header:		
!
!	Content:	
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
!	Author:		Martin Rädel		08.12.2009
!			TU Dresden,Diplomarbeit
!
!	Revision:	
!
! =================================================================================================
!
! Use
!
  use globale_variablen
!
! =================================================================================================
!
  implicit none 
!
! =================================================================================================
!
! Include
!  
#include "finclude/petsc.h"
#include "finclude/petscmat.h"
#include "finclude/petscvec.h"
!
! =================================================================================================
!
! Data types
!
  integer, dimension(num_dof),intent(in) :: num_dof_vec
  integer, intent(in)                    :: dim_dof
  
  Mat, intent(inout)                     :: tMat
  
  integer                                :: node, cid, index, ii, jj, iii, jjj
  PetscInt                               :: pii, pjj
  PetscErrorCode                         :: ierr
  PetscScalar                            :: v
  double precision, dimension(:,:)       :: bnode(3,3)          ! Knotenkoordinaten, die das Koordinatensystem definieren 1 - Ursprung; 2 - z-Achse; 3 - x-z-Ebene
  double precision, dimension(:,:)	 :: BC(3,3)             ! Transformationsmatrix lokal -> global
  
  integer				:: err_code=0
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
  do node = 1, num_nodes
  
    if (nodes(node)%cid .eq. 0) then                            ! keine Transformation
   
      do ii = 1,6                                               ! Hauptdiagonale für diesen Knoten 1-setzen
        
        index = num_dof_vec((node-1)*ndof+ii)
        
        if (index .le. dim_dof) then
          v = 1.d0
          pii = index-1
          call MatSetValue(tMat, pii, pii, v, ADD_VALUES, ierr); CHKERRQ(ierr)
        end if
      
      end do
    
    else
    
      cid          = nodes(node)%cid
      bnode(1,1:3) = nodes(coords(cid)%nids(1))%coords(1:3)
      bnode(2,1:3) = nodes(coords(cid)%nids(2))%coords(1:3)
      bnode(3,1:3) = nodes(coords(cid)%nids(3))%coords(1:3)
      
      call rotation_coordsys (bnode,'lg',BC)
      
      ! Translationen
      do ii = 1,3
        
        iii = num_dof_vec((node-1)*ndof+ii)
        
        if ( iii .LE. dim_dof) then
          do jj = 1,3 
          
            jjj = num_dof_vec((node-1)*ndof+jj)                                         
          
            if (jjj .LE. dim_dof .and. jjj >= iii) then
              pii = ii-1
              pjj = jj-1
              v   = BC(ii,jj)
              call MatSetValue(tMat, pii  , pjj  , v, ADD_VALUES, ierr); CHKERRQ(ierr)  ! Transformation der Translationen
            end if
         
          end do
        end if
      end do
      
      ! Rotationen
      do ii = 1,3
        
        iii = num_dof_vec((node-1)*ndof+ii+3)
        
        if ( iii .LE. dim_dof) then
          do jj = 1,3 
          
            jjj = num_dof_vec((node-1)*ndof+jj+3)                                         
          
            if (jjj .LE. dim_dof) then
              pii = ii-1
              pjj = jj-1
              v   = BC(ii,jj)
              call MatSetValue(tMat, pii  , pjj  , v, ADD_VALUES, ierr); CHKERRQ(ierr)  ! Transformation der Rotationen
            end if
         
          end do
        end if
      end do
    
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
   write(*,*)  			   'transform_nodal_coords'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine transform_nodal_coords
