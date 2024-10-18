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
!> $Id: koordinatensystem_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module koordinatensystem_typen

  implicit none

type coordsyses_type
  type (coord_type), allocatable :: coords (:)
end type coordsyses_type
  
! =================================================================================================
!
! Coordinate System cards
!

type coord_type
  integer                               :: cid      ! Coordinate System-ID           (Integer>0)
  double precision, dimension(3,3)      :: transMat ! Transformationsmatrix lokal zu global. Darin sind die Spalten der Reihe nach der normierte X-, Y- und Z-Vektor (Real)
end type coord_type

logical                                 :: mpi_type_inited = .false.
integer                                 :: mpi_coord_type

contains

! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================

  subroutine ksys_mpi_init()

#include "petsc/finclude/petscsys.h"

#if !defined (PETSC_HAVE_MPIUNI)

    use petscsys
  
    implicit none
    
    INTEGER                             :: count
    INTEGER, dimension(:), allocatable  :: array_of_blockl
    INTEGER, dimension(:), allocatable  :: array_of_types
    INTEGER(kind=MPI_ADDRESS_KIND), dimension(:), allocatable   :: array_of_disp  
    INTEGER(kind=MPI_ADDRESS_KIND), dimension(:), allocatable   :: address
    
    INTEGER                             :: ii
    INTEGER                             :: ierror
    
    TYPE(coord_type)                    :: coord

    !***************************************************************************
    ! coord type
    !***************************************************************************

    count = 2
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! cid
    array_of_types (1) = MPI_INTEGER
    array_of_blockl(1) = 1
    CALL MPI_GET_ADDRESS(coord%cid, address(1), ierror)
    
    ! xAxisVec
    array_of_types (2) = MPI_DOUBLE_PRECISION
    array_of_blockl(2) = 9
    CALL MPI_GET_ADDRESS(coord%transMat, address(2), ierror)
    
    array_of_disp(1) = 0

    DO ii = 2, count

      array_of_disp(ii) = address(ii) - address(1)

    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types,  	   &	
    &				mpi_coord_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_coord_type, ierror)
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)
#endif
  end subroutine ksys_mpi_init
  
! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================  
  
  subroutine bcast_ksys(csyses,is_coord)
#include "petsc/finclude/petscsys.h"
    use petscsys
    use globale_variablen
  
    implicit none
    
    type(coordsyses_type)   :: csyses
    logical                 :: is_coord
    
    PetscMPIInt             :: rank
    PetscErrorCode          :: ierr
    integer                 :: coord_num
    
#if !defined (PETSC_HAVE_MPIUNI)

    if (mpi_type_inited .eq. .false.) then
      call ksys_mpi_init()
      mpi_type_inited = .true.
    end if
    
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
    
    ! coords
    if (is_coord .eq. .true.) then
      if (rank .eq. 0) coord_num = size(csyses%coords,1) 
      call MPI_Bcast (coord_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(csyses%coords(coord_num))
      call MPI_Bcast (csyses%coords, coord_num, MPI_COORD_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if
#endif
  end subroutine bcast_ksys
  
  subroutine free_mem_ksys(koordinatensysteme)
  
    implicit none
    
    type(coordsyses_type) :: koordinatensysteme
    
    if (allocated(koordinatensysteme%coords) .eq. .true.) deallocate(koordinatensysteme%coords)
  
  end subroutine free_mem_ksys

end module koordinatensystem_typen
