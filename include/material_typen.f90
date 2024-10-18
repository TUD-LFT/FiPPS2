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
!> $Id: material_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module material_typen

  implicit none
  
type materials_type
  type (mat1_type), allocatable   :: mat1s (:)
  type (mat8_type), allocatable   :: mat8s (:)
  type (mat20_type), allocatable  :: mat20s (:)
end type materials_type

! =================================================================================================
!
! Material cards
!
! mat1-Remarks:	- mids must be unique for all mat1,...,mat8 entries
! 		- when ym, nu or sm are blank:
! 			- ym and sm my not both be blank
! 			- if ym and nu or sm and nu are both blank, both are set to 0.0
! 			- if only ym or nu or sm are blank the missing item will be computed from
! 			  equation: ym=2*(1+nu)*sm
! 		- if values are specified for all of the properties ym, nu and sm then it is recommended
! 		  that the following relationship is satisfied: abs(1-ym/(2*(1+nu)*G))<0.01
! 		  if this relationship is not desired anisotropic material definition is recommended

type mat1_type
  integer               :: mid  ! Material-ID
  double precision      :: ym   ! Youngs modulus (E-Modul)      (Real>=0.0 or blank)
  double precision      :: sm   ! Shear modulus                 (Real>=0.0 or blank)
  double precision      :: nu   ! Poisson ration                (-1.0<=Real<=0.5 or blank)
  double precision      :: rho  ! Mass density                  (Real>=0.0)
  double precision      :: ath  ! Thermal expansion coefficient (n.i.)
  double precision      :: tref ! Reference temperature for calculation of thermal loads or a temperature-
                                ! dependent ath                 (n.i.)
  double precision      :: ge   ! structural damping coefficient(n.i)
  
  integer, dimension(4) :: fid  ! failure criterion IDs
  integer, dimension(4) :: ifid ! internal failure criterion IDs
  integer, dimension(4) :: iftype ! internal failure criterion Type
end type mat1_type

!--------------------------------------------------------------------------------------------------
!
! Mat8
!
type mat8_type
  integer               :: mid  ! Material-ID
  double precision      :: ym11 ! Youngs modulus (E-Modul)  1-dir   (Real>=0.0 or blank)
  double precision      :: ym22 ! Youngs modulus (E-Modul)  2-dir   (Real>=0.0 or blank)
  double precision      :: nu12 ! Poisson ration            12-dir  (-1.0<=Real<=0.5 or blank)
  double precision      :: sm12 ! Shear modulus             12-dir  (Real>=0.0 or blank)
  double precision      :: sm13 ! Shear modulus             13-dir  (Real>=0.0 or blank)
  double precision      :: sm23 ! Shear modulus             23-dir  (Real>=0.0 or blank)
  double precision      :: rho  ! Mass density                      (n.i.)
  double precision      :: ath11! Thermal expansion coefficient     (n.i.)
  double precision      :: ath22! Thermal expansion coefficient     (n.i.)
  double precision      :: tref ! Reference temperature for calculation of thermal loads or a temperature-
                                ! dependent ath                     (n.i.)
  double precision      :: ge   ! structural damping coefficient    (n.i.)
  
  integer, dimension(4) :: fid  ! failure criterion IDs
  integer, dimension(4) :: ifid ! internal failure criterion IDs
  integer, dimension(4) :: iftype ! internal failure criterion Type
end type mat8_type

!-------------------------------------------------------------------------------
!
! Mat20
!
type mat20_type
  integer               :: mid  ! Material-ID
  double precision      :: ym11 ! Youngs modulus (E-Modul) 1-dir (Real>=0.0 or blank)
  double precision      :: ym22 ! Youngs modulus (E-Modul) 2-dir (Real>=0.0 or blank)
  double precision      :: ym33 ! Youngs modulus (E-Modul) 3-dir (Real>=0.0 or blank) 
  double precision      :: nu12 ! Poisson ration  12-dir  (-1.0<=Real<=0.5 or blank)
  double precision      :: nu13 ! Poisson ration  13-dir  (-1.0<=Real<=0.5 or blank)
  double precision      :: nu23 ! Poisson ration  23-dir  (-1.0<=Real<=0.5 or blank)
  double precision      :: sm12 ! Shear modulus   12-dir  (Real>=0.0 or blank)
  double precision      :: sm13 ! Shear modulus   13-dir  (Real>=0.0 or blank)
  double precision      :: sm23 ! Shear modulus   23-dir  (Real>=0.0 or blank)              
  double precision      :: ath11! Thermal expansion coefficient 1-dir  (Real or blank)
  double precision      :: ath22! Thermal expansion coefficient 2-dir  (Real or blank)
  double precision      :: ath33! Thermal expansion coefficient 3-dir  (Real or blank)
  double precision      :: rho  ! Mass density  (n.i.)
  
  integer, dimension(4) :: fid  ! failure criterion IDs
  integer, dimension(4) :: ifid ! internal failure criterion IDs
  integer, dimension(4) :: iftype ! internal failure criterion Type
end type mat20_type

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
  subroutine bcast_material(materials,is_mat1,is_mat8,is_mat20)
    
#include "petsc/finclude/petscsys.h"
    use petscsys
    use globale_variablen
    
    implicit none
    
    type(materials_type)    :: materials
    logical                 :: is_mat1
    logical                 :: is_mat8
    logical                 :: is_mat20

    PetscMPIInt             :: rank
    PetscErrorCode          :: ierr
    INTEGER                 :: ii
    INTEGER                 :: nummat1
              
#if !defined (PETSC_HAVE_MPIUNI)
              
    call MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr); CHKERRQ(ierr)

    if (is_mat1 .eq. .true.) then
     
      if (rank .eq. 0) nummat1 = size(materials%mat1s,1)
      call MPI_Bcast (nummat1, 1, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      if (rank .ne. 0) allocate(materials%mat1s(nummat1))
      
      do ii = 1,nummat1
      
          call MPI_Bcast (materials%mat1s(ii)%mid,  1,          MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%ym,   1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%sm,   1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%nu,   1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%rho,  1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%ath,  1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%tref, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%ge,   1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
          call MPI_Bcast (materials%mat1s(ii)%fid,  4, MPI_INTEGER, 0, PETSC_COMM_WORLD, ierr); CHKERRQ(ierr)
      
      end do
    
    end if

#endif
    
  end subroutine bcast_material
  
  subroutine free_mem_material(materialien)
  
    implicit none
    
    type(materials_type) :: materialien
    
    if (allocated(materialien%mat1s)  .eq. .true.) deallocate(materialien%mat1s)
    if (allocated(materialien%mat8s)  .eq. .true.) deallocate(materialien%mat8s)
    if (allocated(materialien%mat20s) .eq. .true.) deallocate(materialien%mat20s)
  
  end subroutine free_mem_material

end module material_typen
