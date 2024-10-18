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
!> @brief
!
!
!> @details
!
!> @author
!
!
!> $Id: eigenschaften_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module eigenschaften_typen

  implicit none

type eigenschaften_type
  type (pbeam_type), allocatable    :: pbeams (:)
  type (pshell_type), allocatable   :: pshells (:)
  type (pcomp_type), allocatable    :: pcomps (:)
  type (plsolid_type), allocatable  :: plsolids(:)
end type eigenschaften_type

! =================================================================================================
!
! Property cards
!

!--------------------------------------------------------------------------------------------------
!
! Pbeam
!
! Beam with constant section properties

type pbeam_type
  integer                               :: pid      ! Property-ID                   (Integer>0)
  integer                               :: mid      ! Material-ID                   (Integer>0, must be mat1 type mid)
  double precision                      :: AA       ! Area                          (Real)
  double precision                      :: I11      ! Moment of inertia about z=Izz (Real)
  double precision                      :: I22      ! Moment of inertia about y=Iyy (Real)
  double precision                      :: I12      ! Product of inertia       =Izy (Real, but I11*I22-I12**2>0.0)
  double precision                      :: It       ! Torsional moment of inertia   (Real)
  double precision                      :: t1       ! Thickness in z direction      (Real)
  double precision                      :: t2       ! Thickness in y direction      (Real)
  double precision                      :: angle    ! angle between principal axis  (Real)
                                                    ! and element local coordinate system
!  character(3)                          :: atype    ! angle-type                    ('deg' or 'rad')
  double precision                      :: nsm      ! Nonstruct. mass per unit area (Real)
!internal
  integer                               :: intMat1ID
  double precision                      :: lengthWeight
  double precision                      :: weight
end type pbeam_type

!--------------------------------------------------------------------------------------------------
!
! Pshell

type pshell_type
  integer                               :: pid      ! Property-ID                           (Integer>0)
  integer                               :: mid1     ! Material-ID1 for membrane stiffness   (Integer>0)
  double precision                      :: mt       ! Membrane thickness                    (Real>0)
  integer                               :: mid2     ! Material-ID2 for bending stiffness    (Integer>0)
  double precision                      :: bmr      ! bending moment of inertia ratio 12I/t^3- ratio of the actual bending moment
                                                    ! inertia of the shell I to the bending moment of inertia of a homgeneous shell
                                                    ! t^3/12                                (Real>0.0, Default=1.0)
  integer                               :: mid3     ! Material-ID3 for transverse shear stiffness   (Integer>0)
  double precision                      :: tst      ! transverse shear thickness ratio ts/t - ratio of the shear thickness ts to
                                                    ! the actual (membrane) thickness t - shear correction factor, default value for
                                                    ! homogeneous shell                     (Real>0.0, Default=5.0/6.0)
  double precision                      :: nsm      ! Nonstructural mass per unit area      (Real>0.0)
  double precision                      :: z1       ! fibre distance for stress calculations. the positive direction is determined
                                                    ! by the right-hand rule and the order in which the grid points are listed on the
                                                    ! connection entry                      (Real, Default=-t/2)
  double precision                      :: z2       !                                       (Real, Default= t/2)
  integer                               :: mid4     ! Material-ID4 for membrane-bending-coupling stiffness  (Integer>0)
!internal
  double precision                      :: areaWeight
  double precision                      :: weight
  integer                               :: intMat1ID
  double precision                      :: SOmega
end type pshell_type

!--------------------------------------------------------------------------------------------------
!
! Pcomp

type pcomp_type
  integer                               :: pid      ! Property-ID                   (Integer>0)
  integer                               :: lamid    ! Laminate-ID                   (Integer>0)
  double precision                      :: offset   ! Offset                        (Real)
  double precision                      :: nsm      ! Nonstruct. mass per unit area (Real)
  double precision                      :: sb       ! allowable interlaminar shear stress of the bonding material required if ft
                                                    ! is also specified             (Real>0.0)
  character(2)                          :: ft       ! failure theory                (character or blank)
                                                    ! if blank no failure theory will be performed
!internal (for strain calculation)
  double precision, allocatable, dimension(:,:,:)   :: C        ! Stiffness matrix per layer (unrotated)
  double precision, allocatable, dimension(:,:,:)   :: T        ! Transformation Matrix for C per layer
  double precision, allocatable, dimension(:,:,:)   :: TInv     ! inverse Transformation Matrix for C per layer
  integer, allocatable, dimension(:)                :: nop      ! Number of Integration points over layers (Integer > 0)
  double precision, allocatable, dimension(:)       :: angle    ! Angle of each Layer
  double precision, allocatable, dimension(:)       :: lth      ! Thickness of each layer
  double precision, allocatable, dimension(:,:)     :: ath      ! thermal expansion coefficients for layers in the element coordinate system
  double precision                                  :: tth      ! Total thickness
  integer                                           :: lay      ! number of layers              (Integer>0)
  double precision                                  :: areaWeight
  double precision                                  :: weight
  double precision                                  :: penaltyStiffness ! Auxilary value to calculate the penalty stiffniss for rotation dof 6
  integer, allocatable, dimension(:)                :: intMatID ! internal mat IDs
end type pcomp_type

!--------------------------------------------------------------------------------
!
! Plsolid
!
!

type plsolid_type
  integer                          :: pid           ! Property-ID                   (Integer>0)
  integer                          :: lamid         ! Laminate-ID                   (Integer=>0)
  integer                          :: cid           ! Coordinate System-ID for orientation of
                                                    ! element coordinate system     (Integer=>0, or -1)
  logical                          :: globOut       ! If true, output in global coosy, if false
                                                    ! output in element coosy
  integer                          :: resLay        ! Number of layer for which output strains
                                                    ! and stresses are calculated
!internal
  integer                                         :: lay    ! Number of Layers      (Integer => 0)
  double precision                                :: tth
  double precision, allocatable, dimension(:)     :: lth
  double precision, allocatable, dimension(:)     :: angle
  integer, allocatable, dimension(:)              :: nop    ! Number of Integration points over layers (Integer)
  double precision, allocatable, dimension(:,:,:) :: C
  double precision, allocatable, dimension(:,:)   :: ath
  double precision                                :: weight
  integer, allocatable, dimension(:)              :: intMatID ! internal mat IDs ! internal mat IDs
  
end type plsolid_type

logical                                 :: mpi_type_inited = .false.
integer                                 :: mpi_pbeam_type
integer                                 :: mpi_pshell_type

contains

  subroutine eigensch_mpi_init()

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
    TYPE(pbeam_type)                    :: pbeam
    TYPE(pshell_type)                   :: pshell

    !***************************************************************************
    ! pbeam type
    !***************************************************************************

    count = 12
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! pid
    array_of_types (1) = MPI_INTEGER
    array_of_blockl(1) = 1
    CALL MPI_GET_ADDRESS(pbeam%pid, address(1), ierror)
    
    ! mid
    array_of_types (2) = MPI_INTEGER
    array_of_blockl(2) = 1
    CALL MPI_GET_ADDRESS(pbeam%mid, address(2), ierror)
    
    ! AA
    array_of_types (3) = MPI_DOUBLE_PRECISION
    array_of_blockl(3) = 1
    CALL MPI_GET_ADDRESS(pbeam%AA, address(3), ierror)
    
    ! I11
    array_of_types (4) = MPI_DOUBLE_PRECISION
    array_of_blockl(4) = 1
    CALL MPI_GET_ADDRESS(pbeam%I11, address(4), ierror)
    
    ! I22
    array_of_types (5) = MPI_DOUBLE_PRECISION
    array_of_blockl(5) = 1
    CALL MPI_GET_ADDRESS(pbeam%I22, address(5), ierror)
    
    ! I12
    array_of_types (6) = MPI_INTEGER
    array_of_blockl(6) = 1
    CALL MPI_GET_ADDRESS(pbeam%I12, address(6), ierror)
    
    ! It
    array_of_types (7) = MPI_DOUBLE_PRECISION
    array_of_blockl(7) = 1
    CALL MPI_GET_ADDRESS(pbeam%It, address(7), ierror)
    
    ! t1
    array_of_types (8) = MPI_DOUBLE_PRECISION
    array_of_blockl(8) = 1
    CALL MPI_GET_ADDRESS(pbeam%t1, address(8), ierror)
    
    ! t2
    array_of_types (9) = MPI_DOUBLE_PRECISION
    array_of_blockl(9) = 1
    CALL MPI_GET_ADDRESS(pbeam%t2, address(9), ierror)
    
    ! angle
    array_of_types (10) = MPI_DOUBLE_PRECISION
    array_of_blockl(10) = 1
    CALL MPI_GET_ADDRESS(pbeam%angle, address(10), ierror)
    
    ! nsm
    array_of_types (11) = MPI_DOUBLE_PRECISION
    array_of_blockl(11) = 1
    CALL MPI_GET_ADDRESS(pbeam%nsm, address(11), ierror)
    
    ! intMat1ID
    array_of_types (12) = MPI_INTEGER
    array_of_blockl(12) = 1
    CALL MPI_GET_ADDRESS(pbeam%intMat1ID, address(12), ierror)
    
    array_of_disp(1) = 0

    DO ii = 2, count
       
      array_of_disp(ii) = address(ii) - address(1)

    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types, &
    &                           mpi_pbeam_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_pbeam_type, ierror)
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)
    
    !***************************************************************************
    ! pshell type
    !***************************************************************************

    count = 12
  
    ALLOCATE(array_of_blockl(1:count)); ALLOCATE(array_of_types(1:count))
    ALLOCATE(array_of_disp(1:count)); ALLOCATE(address(1:count))
    
    ! pid
    array_of_types (1) = MPI_INTEGER
    array_of_blockl(1) = 1
    CALL MPI_GET_ADDRESS(pshell%pid, address(1), ierror)
    
    ! mid1
    array_of_types (2) = MPI_INTEGER
    array_of_blockl(2) = 1
    CALL MPI_GET_ADDRESS(pshell%mid1, address(2), ierror)
    
    ! mt
    array_of_types (3) = MPI_DOUBLE_PRECISION
    array_of_blockl(3) = 1
    CALL MPI_GET_ADDRESS(pshell%mt, address(3), ierror)
    
    ! mid2
    array_of_types (4) = MPI_INTEGER
    array_of_blockl(4) = 1
    CALL MPI_GET_ADDRESS(pshell%mid2, address(4), ierror)
    
    ! bmr
    array_of_types (5) = MPI_DOUBLE_PRECISION
    array_of_blockl(5) = 1
    CALL MPI_GET_ADDRESS(pshell%bmr, address(5), ierror)
    
    ! mid3
    array_of_types (6) = MPI_INTEGER
    array_of_blockl(6) = 1
    CALL MPI_GET_ADDRESS(pshell%mid3, address(6), ierror)
    
    ! tst
    array_of_types (7) = MPI_DOUBLE_PRECISION
    array_of_blockl(7) = 1
    CALL MPI_GET_ADDRESS(pshell%tst, address(7), ierror)
    
    ! nsm
    array_of_types (8) = MPI_DOUBLE_PRECISION
    array_of_blockl(8) = 1
    CALL MPI_GET_ADDRESS(pshell%nsm, address(8), ierror)
    
    ! z1
    array_of_types (9) = MPI_DOUBLE_PRECISION
    array_of_blockl(9) = 1
    CALL MPI_GET_ADDRESS(pshell%z1, address(9), ierror)
    
    ! mid4
    array_of_types (10) = MPI_INTEGER
    array_of_blockl(10) = 1
    CALL MPI_GET_ADDRESS(pshell%mid4, address(10), ierror)
    
    ! intMat1ID
    array_of_types (11) = MPI_INTEGER
    array_of_blockl(11) = 1
    CALL MPI_GET_ADDRESS(pshell%intMat1ID, address(11), ierror)
    
    ! SOmega
    array_of_types (12) = MPI_DOUBLE_PRECISION
    array_of_blockl(12) = 1
    CALL MPI_GET_ADDRESS(pshell%SOmega, address(12), ierror)
    
    array_of_disp(1) = 0

    DO ii = 2, count
       
      array_of_disp(ii) = address(ii) - address(1)

    END DO

    CALL MPI_TYPE_CREATE_STRUCT(count, array_of_blockl, array_of_disp, array_of_types,  	   &	
    &				mpi_pshell_type, ierror)
    
    CALL MPI_TYPE_COMMIT(mpi_pshell_type, ierror)
    
    DEALLOCATE(array_of_blockl, array_of_types, array_of_disp, address)

#endif

  end subroutine eigensch_mpi_init

  subroutine bcast_eigenschaften(eigenschaften,is_pbeam,is_pshell,is_pcomp,is_plsolid)
  
#include "petsc/finclude/petscsys.h"
    use petscsys
  
    implicit none
    
    type(eigenschaften_type) :: eigenschaften
    logical                  :: is_pbeam
    logical                  :: is_pshell
    logical                  :: is_pcomp,is_plsolid
    
    PetscMPIInt    :: rank
    PetscErrorCode :: ierr
    integer        :: pbeam_num, pshell_num, pcomp_num, plsolid_num

    !
    ! zum Übertragen der plsolids
    !
    integer                                      :: ii, jj, mm, nn, pos
    integer                                      :: nlay
    double precision, allocatable, dimension(:)  :: C_vec, T_vec, TInv_vec, ath_vec
    
#if !defined (PETSC_HAVE_MPIUNI)
    
    if (mpi_type_inited .eq. .false.) then
      call eigensch_mpi_init()
      mpi_type_inited = .true.
    end if
                    
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)
    
    ! pbeams
    if (is_pbeam .eq. .true.) then
      if (rank .eq. 0) pbeam_num = size(eigenschaften%pbeams,1) 
      call MPI_Bcast (pbeam_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(eigenschaften%pbeams(pbeam_num))
      call MPI_Bcast (eigenschaften%pbeams, pbeam_num, MPI_PBEAM_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if
    
    ! pshells
    if (is_pshell .eq. .true.) then
      if (rank .eq. 0) pshell_num = size(eigenschaften%pshells,1) 
      call MPI_Bcast (pshell_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(eigenschaften%pshells(pshell_num))
      call MPI_Bcast (eigenschaften%pshells, pshell_num, MPI_PSHELL_TYPE, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
    end if
    
    ! pcomps
    if (is_pcomp .eq. .true.) then
      ! Hier kann nicht über einen MPI-Type übertragen werden, da jedes allokierbare Array eine separate Länge in Abhängigkeit von der Lagenanzahl hat
      if (rank .eq. 0) pcomp_num = size(eigenschaften%pcomps,1) 
      call MPI_Bcast (     pcomp_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(eigenschaften%pcomps(pcomp_num))
      
      do ii = 1, pcomp_num
        if (rank .EQ. 0) nlay = eigenschaften%pcomps(ii)%lay
        call MPI_Bcast(nlay, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
        
        allocate(C_vec(36*nlay))
        allocate(T_vec(36*nlay))
        allocate(TInv_vec(36*nlay))
        allocate(ath_vec(6*nlay))
      
        if (rank .EQ. 0) then
          C_vec = reshape( eigenschaften%pcomps(ii)%C, (/ 36*nlay /) )
          T_vec = reshape( eigenschaften%pcomps(ii)%T, (/ 36*nlay /) )
          TInv_vec = reshape( eigenschaften%pcomps(ii)%TInv, (/ 36*nlay /) )
          ath_vec = reshape( eigenschaften%pcomps(ii)%ath, (/ 6*nlay /) )
        else
          allocate(eigenschaften%pcomps(ii)%lth(nlay))
          allocate(eigenschaften%pcomps(ii)%nop(nlay))
          allocate(eigenschaften%pcomps(ii)%angle(nlay))
          allocate(eigenschaften%pcomps(ii)%intMatID(nlay))
        end if
   
        call MPI_Bcast(eigenschaften%pcomps(ii)%pid,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%lamid,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%offset,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%nsm,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%sb,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%ft,2,MPI_CHARACTER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        
        call MPI_Bcast(C_vec,36*nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(T_vec,36*nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(TInv_vec,36*nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%nop,nlay,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%angle,nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%lth,nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(ath_vec,6*nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%tth,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%lay,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        
        call MPI_Bcast(eigenschaften%pcomps(ii)%areaWeight,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%weight,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%penaltyStiffness,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%pcomps(ii)%intMatID,nlay,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
      
        if (rank .NE. 0) then
          allocate(eigenschaften%pcomps(ii)%C(6,6,nlay))
          eigenschaften%pcomps(ii)%C = reshape( C_vec, (/ 6,6,nlay /) )
          allocate(eigenschaften%pcomps(ii)%T(6,6,nlay))
          eigenschaften%pcomps(ii)%T = reshape( T_vec, (/ 6,6,nlay /) )
          allocate(eigenschaften%pcomps(ii)%TInv(6,6,nlay))
          eigenschaften%pcomps(ii)%TInv = reshape( TInv_vec, (/ 6,6,nlay /) )
          allocate(eigenschaften%pcomps(ii)%ath(3,nlay))
          eigenschaften%pcomps(ii)%ath = reshape( ath_vec, (/ 6,nlay /) )
        end if
     
        deallocate(C_vec, T_vec, TInv_vec, ath_vec)
        
      end do
    end if   
    
    if (is_plsolid .EQ. .TRUE.) then 
      ! Hier kann nicht über einen MPI-Type übertragen werden, da jedes allokierbare Array eine separate Länge in Abhängigkeit von der Lagenanzahl hat
      
      ! Die Felder werden hier vorher deallokiert, da im Rahmen einer Multisteprechnung eine Aktualisierung der Eigenschaften notwendig sein kann
      
      if (rank .EQ. 0) plsolid_num = size(eigenschaften%plsolids,1)
      call MPI_Bcast(plsolid_num, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .NE. 0) then
        if (allocated(eigenschaften%plsolids) .eq. .true.) deallocate(eigenschaften%plsolids)
        allocate(eigenschaften%plsolids(plsolid_num))
      end if
      
      do ii = 1, plsolid_num
        if (rank .EQ. 0) nlay = eigenschaften%plsolids(ii)%lay
        call MPI_Bcast(nlay, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      
        allocate(C_vec(36*nlay))
        allocate(ath_vec(6*nlay))
      
        if (rank .EQ. 0) then
          do jj = 1, nlay
            do mm = 1, 6
              pos = (jj-1)*6+mm
              ath_vec(pos) = eigenschaften%plsolids(ii)%ath(mm,jj)
              do nn = 1, 6
                pos = (jj-1)*36+(mm-1)*6+nn
                C_vec(pos) = eigenschaften%plsolids(ii)%C(nn,mm,jj)
              end do
            end do
          end do
        else
          if (allocated(eigenschaften%plsolids(ii)%lth) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%lth)
          if (allocated(eigenschaften%plsolids(ii)%nop) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%nop)
          if (allocated(eigenschaften%plsolids(ii)%angle) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%angle)
          if (allocated(eigenschaften%plsolids(ii)%intMatID) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%intMatID)
          allocate(eigenschaften%plsolids(ii)%lth(nlay))
          allocate(eigenschaften%plsolids(ii)%nop(nlay))
          allocate(eigenschaften%plsolids(ii)%angle(nlay))
          allocate(eigenschaften%plsolids(ii)%intMatID(nlay))
        end if
   
        call MPI_Bcast(eigenschaften%plsolids(ii)%pid,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%lamid,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%cid,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%globOut,1,MPI_LOGICAL,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%resLay,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%lay,1,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%tth,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%lth,nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%nop,nlay,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%angle,nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(C_vec,36*nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(ath_vec,6*nlay,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%weight,1,MPI_DOUBLE_PRECISION,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
        call MPI_Bcast(eigenschaften%plsolids(ii)%intMatID,nlay,MPI_INTEGER,0,PETSC_COMM_WORLD,ierr); CHKERRQ(ierr)
      
        if (rank .NE. 0) then
          if (allocated(eigenschaften%plsolids(ii)%ath) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%ath)
          if (allocated(eigenschaften%plsolids(ii)%C) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%C)
          allocate(eigenschaften%plsolids(ii)%ath(6,nlay))
          allocate(eigenschaften%plsolids(ii)%C(6,6,nlay))
      
          do jj = 1, nlay
            do mm = 1, 6
              pos = (jj-1)*6+mm
              eigenschaften%plsolids(ii)%ath(mm,jj) = ath_vec(pos)
              do nn = 1, 6
                pos = (jj-1)*36+(mm-1)*6+nn
                eigenschaften%plsolids(ii)%C(nn,mm,jj) = C_vec(pos)
              end do
            end do
          end do
        end if
     
        deallocate(C_vec, ath_vec)
     
      end do
    end if
#endif
  end subroutine bcast_eigenschaften
  
  subroutine free_mem_eigenschaften(eigenschaften)
  
    implicit none
    
    type(eigenschaften_type) :: eigenschaften
    integer :: ii
    
    if (allocated(eigenschaften%pshells)   .eq. .true.) deallocate(eigenschaften%pshells)
    if (allocated(eigenschaften%pbeams)    .eq. .true.) deallocate(eigenschaften%pbeams)
    if (allocated(eigenschaften%pcomps)    .eq. .true.) then
      do ii = 1, size(eigenschaften%pcomps,1)
        if (allocated(eigenschaften%pcomps(ii)%C) .eq. .true.)        deallocate(eigenschaften%pcomps(ii)%C)
        if (allocated(eigenschaften%pcomps(ii)%T) .eq. .true.)        deallocate(eigenschaften%pcomps(ii)%T)
        if (allocated(eigenschaften%pcomps(ii)%TInv) .eq. .true.)     deallocate(eigenschaften%pcomps(ii)%TInv)
        if (allocated(eigenschaften%pcomps(ii)%nop) .eq. .true.)      deallocate(eigenschaften%pcomps(ii)%nop)
        if (allocated(eigenschaften%pcomps(ii)%angle) .eq. .true.)    deallocate(eigenschaften%pcomps(ii)%angle)
        if (allocated(eigenschaften%pcomps(ii)%lth) .eq. .true.)      deallocate(eigenschaften%pcomps(ii)%lth)
        if (allocated(eigenschaften%pcomps(ii)%ath) .eq. .true.)      deallocate(eigenschaften%pcomps(ii)%ath)
        if (allocated(eigenschaften%pcomps(ii)%intMatID) .eq. .true.) deallocate(eigenschaften%pcomps(ii)%intMatID)
      end do
      deallocate(eigenschaften%pcomps)
    end if
    if (allocated(eigenschaften%plsolids)  .eq. .true.) then
      do ii = 1, size(eigenschaften%plsolids,1)
        if (allocated(eigenschaften%plsolids(ii)%lth) .eq. .true.)   deallocate(eigenschaften%plsolids(ii)%lth)
        if (allocated(eigenschaften%plsolids(ii)%angle) .eq. .true.) deallocate(eigenschaften%plsolids(ii)%angle)
        if (allocated(eigenschaften%plsolids(ii)%nop) .eq. .true.)   deallocate(eigenschaften%plsolids(ii)%nop)
        if (allocated(eigenschaften%plsolids(ii)%C) .eq. .true.)     deallocate(eigenschaften%plsolids(ii)%C)
        if (allocated(eigenschaften%plsolids(ii)%ath) .eq. .true.)   deallocate(eigenschaften%plsolids(ii)%ath)
      end do
      deallocate(eigenschaften%plsolids)
    end if
  
  end subroutine free_mem_eigenschaften
  
end module eigenschaften_typen
