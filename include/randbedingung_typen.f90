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
!> $Id: randbedingung_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module randbedingung_typen

  implicit none
  
type randbedinungs_type
  type (spcadd_type), allocatable :: spcadds (:)
  type (mpcadd_type), allocatable :: mpcadds (:)
  type (spc1_type), allocatable :: spc1s (:)
  type (spcd_type), allocatable :: spcds (:)
  type (mpc_type), allocatable :: mpcs (:)
  type (couplings_type), allocatable :: couplings (:)
  type (contact_node_beam2_type), allocatable :: contact_node_beam2(:)
  type (contact_node_quad8_type), allocatable :: contact_node_quad8(:)
  type (contact_node_lsolid20_type), allocatable :: contact_node_lsolid20(:)
end type randbedinungs_type

! =================================================================================================
!
! Constraint Cards
!

! Constraint sets

type spcadd_type
  integer               :: scid     ! Spc-set-ID                    (Integer>0)
  integer               :: sid      ! Spc-ID                        (Integer>0)
end type spcadd_type

type mpcadd_type
  integer               :: scid     ! Mpc-set-ID                    (Integer>0)
  integer               :: sid      ! Mpc-ID                        (Integer>0)
end type mpcadd_type

! Single-point constraint

type spc1_type
  integer               :: sid      ! Spc-ID                        (Integer>0)
  integer               :: dof      ! Degrees of freedom to be locked (Integer>0)
  integer               :: n1       ! Node                          (Integer>0)
  logical               :: thru     ! Multiple nodes with this spc  (.true. or .false.)
  integer               :: nn       ! End node with this spc        (Integer>0)
end type spc1_type

! Enforced motion value

type spcd_type

   integer               :: sid     ! Spcd-ID                       (Integer>0)
   integer               :: nid     ! Node-ID                       (Integer>0)
   integer               :: dof     ! Degree of freedom spcd is applied on (1<=Integer<=6)
   double precision      :: val     ! value of enforced motion      (Real)
   
end type

! Multi-point constraint

! Auxilary Type for MPCs
type mpc_dof_type
  integer               :: nid      ! Node-ID                       (Integer>0)
  integer               :: dof      ! Degree of Freedom             (Integer>0)
  double precision      :: fac      ! Dactor for the DOF            (Real)
end type mpc_dof_type

type mpc_type
  integer               :: sid      ! MPC-ID                        (Integer>0)
  integer               :: mpc_type	! MPC-Type                      (0 - external, 1 - internal)
  type(mpc_dof_type)    :: dependend	! Dependend DOF            (mpc_dof_type)
  type(mpc_dof_type), allocatable :: independend(:)! Independend DOFs (mpc_dof_type(:))
end type mpc_type

! =================================================================================================
!
! Couplings cards
!

type couplings_type
  integer               :: cpsid    ! CouplingSet-ID                (Integer>0)
  integer               :: dof      ! Degree of Freedom             (Integer>0)
  integer               :: nid      ! Node-ID                       (Integer>0)
end type couplings_type

type contact_node_beam2_type
  integer                        :: beam2ID
  integer, allocatable           :: nodeIDs(:)   ! Node-IDs to be connected with this element
  double precision, allocatable  :: xi(:)        ! Elementkoordinate
end type contact_node_beam2_type

type contact_node_quad8_type
  integer               :: quad8ID
  integer, allocatable  :: nodeIDs(:)   ! Node-IDs to be connected with this element
end type contact_node_quad8_type

type contact_node_lsolid20_type
  integer               :: lsolid20ID
  integer, allocatable  :: nodeIDs(:)   ! Node-IDs to be connected with this element
end type contact_node_lsolid20_type

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
  subroutine bcast_rb(randbedingungen,is_mpc)
    
#include "petsc/finclude/petscsys.h"
    use petscsys
    use globale_variablen
    
    implicit none
    
    type(randbedinungs_type)    :: randbedingungen
    logical                     :: is_mpc
    
    PetscMPIInt                 :: rank
    PetscErrorCode              :: ierr
    INTEGER                     :: ii, jj
    INTEGER                     :: nummpc, numindependend

#if !defined (PETSC_HAVE_MPIUNI)
                    
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

    if (is_mpc .eq. .true.) then
     
      if (rank .eq. 0) nummpc = size(randbedingungen%mpcs,1)
      call MPI_Bcast (nummpc, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(randbedingungen%mpcs(nummpc))
      
      do ii = 1,nummpc
      
        call MPI_Bcast (randbedingungen%mpcs(ii)%sid          , 1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (randbedingungen%mpcs(ii)%mpc_type     , 1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (randbedingungen%mpcs(ii)%dependend%nid, 1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (randbedingungen%mpcs(ii)%dependend%dof, 1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        call MPI_Bcast (randbedingungen%mpcs(ii)%dependend%fac, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)

        if (rank .eq. 0) numindependend = size(randbedingungen%mpcs(ii)%independend,1)
        call MPI_Bcast (numindependend, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        if (rank .ne. 0) allocate(randbedingungen%mpcs(ii)%independend(numindependend))

        do jj = 1,numindependend
          call MPI_Bcast (randbedingungen%mpcs(ii)%independend(jj)%nid, 1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (randbedingungen%mpcs(ii)%independend(jj)%dof, 1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (randbedingungen%mpcs(ii)%independend(jj)%fac, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        end do
      
      end do
    
    end if

#endif

  end subroutine bcast_rb

! =================================================================================================
!
!> @brief
!
!> @details
!
!> @author 
!
! =================================================================================================  
  subroutine free_mem_rb(randbedingungen)
  
    implicit none
    
    type(randbedinungs_type) :: randbedingungen
    integer :: ii
    
    if (allocated(randbedingungen%spcadds)   .eq. .true.) deallocate(randbedingungen%spcadds)
    if (allocated(randbedingungen%mpcadds)   .eq. .true.) deallocate(randbedingungen%mpcadds)
    if (allocated(randbedingungen%spc1s)     .eq. .true.) deallocate(randbedingungen%spc1s)
    if (allocated(randbedingungen%spcds)     .eq. .true.) deallocate(randbedingungen%spcds)
    if (allocated(randbedingungen%mpcs)      .eq. .true.) then
      do ii = 1, size(randbedingungen%mpcs,1)
        deallocate(randbedingungen%mpcs(ii)%independend)
      end do
      deallocate(randbedingungen%mpcs)
    end if
    if (allocated(randbedingungen%couplings) .eq. .true.) deallocate(randbedingungen%couplings)
    if (allocated(randbedingungen%contact_node_beam2)  .eq. .true.) then
      do ii = 1,size(randbedingungen%contact_node_beam2,1)
        deallocate(randbedingungen%contact_node_beam2(ii)%nodeIDs)
        deallocate(randbedingungen%contact_node_beam2(ii)%xi)
      end do
      deallocate(randbedingungen%contact_node_beam2)
    end if
    if (allocated(randbedingungen%contact_node_quad8)  .eq. .true.) then
      do ii = 1,size(randbedingungen%contact_node_quad8,1)
        deallocate(randbedingungen%contact_node_quad8(ii)%nodeIDs)
      end do
      deallocate(randbedingungen%contact_node_quad8)
    end if
    if (allocated(randbedingungen%contact_node_lsolid20)  .eq. .true.) then
      do ii = 1,size(randbedingungen%contact_node_lsolid20,1)
        deallocate(randbedingungen%contact_node_lsolid20(ii)%nodeIDs)
      end do
      deallocate(randbedingungen%contact_node_lsolid20)
    end if
  
  end subroutine free_mem_rb

end module randbedingung_typen
