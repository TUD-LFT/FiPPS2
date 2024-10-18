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
!> $Id: last_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module last_typen

  implicit none
  
type lasten_type
  type (load_type), allocatable         :: loads (:)
  type (force_type), allocatable        :: forces (:)
  type (moment_type), allocatable       :: moments (:)
  type (p2load_type), allocatable       :: p2loads (:)
  type (p8load_type), allocatable       :: p8loads (:)
  type (p20load_type), allocatable      :: p20loads (:)
  type (aeroload2d_type), allocatable   :: aeroload2ds (:)
  type (aeroload3d_type), allocatable   :: aeroload3ds (:)
  type (nodetemp_type), allocatable     :: nodetemps (:)
  type (elementtemp_type), allocatable  :: beam2temps (:)
  type (elementtemp_type), allocatable  :: quad8temps (:)
  type (elementtemp_type), allocatable  :: lsolid20temps (:)
  type (subcase_type), allocatable      :: subcases (:)
end type lasten_type

!
! =================================================================================================
!
! Load cards

type load_type
  integer               :: lcid     ! Load-case-ID                  (Integer>0)
  double precision      :: sfac     ! Scaling factor for all loads  (Real)
  double precision      :: sfaci    ! Scaling factor for load 1     (Real)
  integer               :: lidi     ! Load-ID 1                     (Integer>0)
end type load_type

!--------------------------------------------------------------------------------------------------
!
! Force
!
type force_type
  integer               :: lid      ! Load-set-ID                   (Integer>0)
  integer               :: nid      ! Node-ID                       (Integer>0)
  double precision      :: fac      ! Scaling factor                (Real)
  double precision, dimension(1:3)  :: ni   ! Components of vector force is acting, measured in coordinate system
                                            ! defined by CID        (Real, at last 1 Ni /= 0)
end type force_type

!--------------------------------------------------------------------------------------------------
!
! Moment
!
type moment_type
  integer               :: lid      ! Load-set-ID                   (Integer>0)
  integer               :: nid      ! Node-ID                       (Integer>0)
  double precision      :: fac      ! Scaling factor                (Real)
  double precision, dimension(1:3)  :: ni   ! Components of vector moment is acting, measured in coordinate system
                                            ! defined by CID        (Real, at last 1 Ni /= 0)
end type moment_type

!--------------------------------------------------------------------------------------------------
!
! P2load
!
type p2load_type                            ! constant line load on beam2-type-Element
  integer                           :: lid  ! Load-set-ID           (Integer>0)
  integer                           :: eid1 ! Element-ID            (Integer>0)
  integer                           :: dir  ! ID of the direction in which the line load is applied
  double precision, dimension(1:2)  :: pi   ! Line load at nodes 1-2, only constant line load P1=P2 implemented
                                            ! right now, defined at P1	(Real, at last 1 Pi /= 0)
  logical                           :: thru ! true if multiple elements
  integer                           :: eid2 ! end of elements with pload	(Integer, eid1<eid2)
end type p2load_type
!--------------------------------------------------------------------------------------------------
!
! p8load
!
type p8load_type                ! constant pressure on quad4-type-Element
  integer               :: lid  ! Load-set-ID                       (Integer>0)
  integer               :: eid1 ! Element-ID                        (Integer>0)
  integer               :: cid  ! Coordinate-System-ID              (Integer>=0, Default 0)
  double precision, dimension(1:4)  :: pi   ! pressure at corner nodes 1-4, only constant pressure P1=P2=P3=P4 implemented
                                            ! right now, defined at P1	(Real, at last 1 Pi /= 0)
  logical               :: thru ! true if multiple elements
  integer               :: eid2 ! end of elements with pload        (Integer, eid1<eid2)
end type p8load_type

!--------------------------------------------------------------------------------------------------
!
! p20load
!
type p20load_type               ! constant pressure on solid20-type-Element
  integer               :: lid  ! Load-set-ID           (Integer>0)
  integer               :: eid1 ! Element-ID            (Integer>0)
  double precision      :: p    ! pressure value
  integer               :: surf ! Surface-ID where pressure is applied (1,2,...,6 allowed)
  logical               :: thru ! true if multiple elements
  integer               :: eid2 ! end of elements with pload    (Integer, eid1<eid2)
end type p20load_type

!--------------------------------------------------------------------------------------------------
!
! aeroload2d
!
type aeroload2d_type            ! aerodynamic load from 2d panel methods
  integer               :: lid  ! Load-set-ID           (Integer>0)
  integer               :: mthd ! Aerodynamic mehtod    (Integer>0)
  double precision      :: dfac ! damping factor for displacements during fluid structure coupling
                                !                       (Real > 0.0)
end type aeroload2d_type

!--------------------------------------------------------------------------------------------------
!
! aeroload3d
!
type aeroload3d_type            ! aerodynamic load from 3d panel methods
  integer               :: lid  ! Load-set-ID           (Integer>0)
  integer               :: mthd ! Aerodynamic mehtod    (Integer>0)
  double precision      :: dfac ! damping factor for displacements during fluid structure coupling
                                !                       (Real > 0.0)
end type aeroload3d_type

! Temperature at node

type nodetemp_type
  integer                :: lid  ! Load-set-ID                           (Integer>0)
  integer                :: nid  ! Node-ID                               (Integer>0)
  double precision       :: temp ! Temperature (T-T0)                    (Real)
end type nodetemp_type

! Temperature at elemente

type elementtemp_type
  integer                :: lid  ! Load-set-ID                           (Integer>0)
  integer                :: eid  ! Element-ID                            (Integer>0)
  double precision       :: temp ! Temperature (T-T0)                    (Real)
end type elementtemp_type

! =================================================================================================
!
! Subcase cards
!

type subcase_type
  integer               :: scid     ! Subcase-ID                         (Integer>0)
  integer               :: spcaddid ! Spcadd-ID                          (Integer>0)
  integer               :: loadid   ! Load-ID                            (Integer>0)
  integer               :: mpcaddid ! Mpcadd-ID                          (Integer>0)
  logical               :: skipBuckling ! Skip Buckling Flag             (Logical .TRUE. or .FALSE.)
  logical               :: upgeom   ! Update geometry                    (Logical .TRUE. or .FALSE.)
  logical               :: upstress ! Update stresses                    (Logical .TRUE. or .FALSE.)
  logical               :: output   ! Update stresses                    (Logical .TRUE. or .FALSE.)
  logical               :: readApameInput
  
  !Internal
  logical               :: upmats   ! Update matricies                   (Logical .TRUE. or .FALSE.)
  integer               :: aeroloadID  ! Internal ID for aeroload
  double precision      :: aeroloadFac ! Load factor for aeroload
end type subcase_type

!
! =================================================================================================

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

  subroutine bcast_lasten(lasten)
    
#include "petsc/finclude/petscsys.h"
    use petscsys
    
    implicit none

    PetscMPIInt     :: rank
    PetscErrorCode  :: ierr
    type(lasten_type) :: lasten
    INTEGER         :: ii
    INTEGER         :: numsub
    
#if !defined (PETSC_HAVE_MPIUNI)
    
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
    
    if (rank .eq. 0) numsub = size(lasten%subcases,1)
    call MPI_Bcast (numsub, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    if (rank .ne. 0) allocate(lasten%subcases(numsub))
      
    do ii = 1,numsub
    
        call MPI_Bcast (lasten%subcases(ii)%scid,          1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%spcaddid,      1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%loadid,        1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%mpcaddid,      1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%skipBuckling,  1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%upgeom,        1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%upstress,      1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%output,        1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (lasten%subcases(ii)%upmats,        1, MPI_LOGICAL, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    
    end do
    
#endif
    
  end subroutine bcast_lasten
  
  subroutine free_mem_lasten(lasten)
  
    implicit none
    
    type(lasten_type) :: lasten
    
    if (allocated(lasten%loads)         .eq. .true.) deallocate(lasten%loads)
    if (allocated(lasten%forces)        .eq. .true.) deallocate(lasten%forces)
    if (allocated(lasten%moments)       .eq. .true.) deallocate(lasten%moments)
    if (allocated(lasten%p2loads)       .eq. .true.) deallocate(lasten%p2loads)
    if (allocated(lasten%p8loads)       .eq. .true.) deallocate(lasten%p8loads)
    if (allocated(lasten%p20loads)      .eq. .true.) deallocate(lasten%p20loads)
    if (allocated(lasten%aeroload2ds)   .eq. .true.) deallocate(lasten%aeroload2ds)
    if (allocated(lasten%aeroload3ds)   .eq. .true.) deallocate(lasten%aeroload3ds)
    if (allocated(lasten%nodetemps)     .eq. .true.) deallocate(lasten%nodetemps)
    if (allocated(lasten%beam2temps)    .eq. .true.) deallocate(lasten%beam2temps)
    if (allocated(lasten%quad8temps)    .eq. .true.) deallocate(lasten%quad8temps)
    if (allocated(lasten%lsolid20temps) .eq. .true.) deallocate(lasten%lsolid20temps)
    if (allocated(lasten%subcases)      .eq. .true.) deallocate(lasten%subcases)
  
  end subroutine free_mem_lasten
  
end module last_typen
