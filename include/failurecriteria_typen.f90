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
!> @details
!
!> @author
!
!> $Id: failurecriteria_typen.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
module failurecriteria_typen

  implicit none
  
  
type failurecriteria_type
  type (fail_tresca_type), allocatable          :: failTrescas (:)
  type (fail_mises_type), allocatable           :: failMises (:)
  type (fail_maxprincstress_type), allocatable  :: failMaxprincstresses (:)
  type (fail_puck_type), allocatable            :: failPucks (:)
  type (fail_hill_type), allocatable            :: failHills (:)
  type (fail_norris_type), allocatable          :: failNorris (:)
  type (fail_fibre_type), allocatable           :: failFibres (:)
  type (fail_maxstrain_type), allocatable       :: failMaxStrains (:)
  type (fail_cuntze_type), allocatable          :: failCuntzes (:)
  type (fail_maxstrain3d_type), allocatable     :: failMaxStrain3Ds (:)
  type (fail_tsaiwu3d_type), allocatable        :: failTsaiWu3Ds (:)
end type failurecriteria_type

! =================================================================================================
!
! Failure criterion cards
!

!--------------------------------------------------------------------------------------------------
!
! 2D Tresca failure criterion
!
type fail_tresca_type
  integer               :: fid  ! failure criterion ID
  double precision      :: ys   ! yield stress                  (Real>=0.0)
end type fail_tresca_type

!--------------------------------------------------------------------------------------------------
!
! 2D von Mises failure criterion
!
type fail_mises_type
  integer               :: fid  ! failure criterion ID
  double precision      :: ys   ! yield stress                  (Real>=0.0)
end type fail_mises_type

!--------------------------------------------------------------------------------------------------
!
! 2D maximum principal Stresses failure criterion
!
type fail_maxprincstress_type
  integer               :: fid  ! failure criterion ID
  double precision      :: ys   ! yield stress                  (Real>=0.0)
  double precision      :: ysC  ! yield stress compression      (Real>=0.0)
end type fail_maxprincstress_type

!--------------------------------------------------------------------------------------------------
!
! 2D Puck failure criterion
!
type fail_puck_type
  integer               :: fid  ! failure criterion ID
  double precision      :: RParTen ! Yield Stress in fibre direction tension           (Real>=0.0 or blank)
  double precision      :: RParCom ! Yield Stress in fibre direction compression       (Real>=0.0 or blank)
  double precision      :: RNorTen ! Yield Stress normal to fibre direction tension    (Real>=0.0 or blank)
  double precision      :: RNorCom ! Yield Stress normal to fibre direction compression(Real>=0.0 or blank)
  double precision      :: RShear  ! Yield Stress shear                                (Real>=0.0 or blank)
  double precision      :: Pspd    
  double precision      :: Pspz
  double precision      :: a0
  double precision      :: lambdamin
end type fail_puck_type

!--------------------------------------------------------------------------------------------------
!
! 2D Hill failure criterion
!
type fail_hill_type
  integer               :: fid  ! failure criterion ID
  double precision      :: RParTen ! Yield Stress in fibre direction tension           (Real>=0.0 or blank)
  double precision      :: RParCom ! Yield Stress in fibre direction compression       (Real>=0.0 or blank)
  double precision      :: RNorTen ! Yield Stress normal to fibre direction tension    (Real>=0.0 or blank)
  double precision      :: RNorCom ! Yield Stress normal to fibre direction compression(Real>=0.0 or blank)
  double precision      :: RShear  ! Yield Stress shear                                (Real>=0.0 or blank)
  double precision      :: F12star
end type fail_hill_type

!--------------------------------------------------------------------------------------------------
!
! 2D Norris failure criterion
!
type fail_norris_type
  integer               :: fid  ! failure criterion ID
  double precision      :: RPar    ! Failure Stress in fibre direction (blunt notch)          (Real>=0.0 or blank)
  double precision      :: RNor    ! Failure Stress normal to fibre direction (blunt notch)   (Real>=0.0 or blank)
  double precision      :: RShear  ! Failure Stress shear (blunt notch)                       (Real>=0.0 or blank)
end type fail_norris_type

!--------------------------------------------------------------------------------------------------
!
! Fibre failure criterion
!
type fail_fibre_type
  integer               :: fid  ! failure criterion ID
  double precision      :: RParTen ! Yield Stress in fibre direction tension    (Real>=0.0 or blank)
  double precision      :: RParCom ! Yield Stress in fibre direction compression(Real>=0.0 or blank)
end type fail_fibre_type

!--------------------------------------------------------------------------------------------------
!
! Maximum strain failure criterion
!
type fail_maxstrain_type
  integer               :: fid  ! failure criterion ID
  double precision      :: epsParTen ! Failure strain in fibre direction tension           (Real>=0.0 or blank)
  double precision      :: epsParCom ! Failure strain in fibre direction compression       (Real>=0.0 or blank)
  double precision      :: epsNorTen ! Failure strain normal to fibre direction tension    (Real>=0.0 or blank)
  double precision      :: epsNorCom ! Failure strain normal to fibre direction compression(Real>=0.0 or blank)
  double precision      :: epsShear  ! Failure strain shear                                (Real>=0.0 or blank)
  logical               :: useGlobal ! use strains in global coordinate system             (Logical)
end type fail_maxstrain_type

!--------------------------------------------------------------------------------------------------
!
! 2D Cuntze failure criterion
!
type fail_cuntze_type
  integer               :: fid  ! failure criterion ID
  double precision      :: RParTen ! Yield Stress in fibre direction tension           (Real>=0.0 or blank)
  double precision      :: RParCom ! Yield Stress in fibre direction compression       (Real>=0.0 or blank)
  double precision      :: RNorTen ! Yield Stress normal to fibre direction tension    (Real>=0.0 or blank)
  double precision      :: RNorCom ! Yield Stress normal to fibre direction compression(Real>=0.0 or blank)
  double precision      :: RShear  ! Yield Stress shear                                (Real>=0.0 or blank)
  double precision      :: muNorPar
  double precision      :: m
end type fail_cuntze_type

!--------------------------------------------------------------------------------------------------
!
! 3D Maximum strain failure criterion
!
type fail_maxstrain3d_type
  integer               :: fid  ! failure criterion ID
  double precision      :: eps11Ten   ! Failure strain in 11-direction tension              (Real>=0.0 or blank)
  double precision      :: eps11Com   ! Failure strain in 11-direction compression          (Real>=0.0 or blank)
  double precision      :: eps22Ten   ! Failure strain in 22-direction tension              (Real>=0.0 or blank)
  double precision      :: eps22Com   ! Failure strain in 22-direction compression          (Real>=0.0 or blank)
  double precision      :: eps33Ten   ! Failure strain in 33-direction tension              (Real>=0.0 or blank)
  double precision      :: eps33Com   ! Failure strain in 33-direction compression          (Real>=0.0 or blank)
  double precision      :: eps12Shear ! Failure strain in 12-direction shear                (Real>=0.0 or blank)
  double precision      :: eps13Shear ! Failure strain in 13-direction shear                (Real>=0.0 or blank)
  double precision      :: eps23Shear ! Failure strain in 23-direction shear                (Real>=0.0 or blank)
  logical               :: useGlobal  ! use strains in element coordinate system            (Logical)
end type fail_maxstrain3d_type

!--------------------------------------------------------------------------------------------------
!
! 3D Maximum strain failure criterion
!
type fail_tsaiWu3d_type
  integer               :: fid  ! failure criterion ID
  double precision      :: R11Ten   ! Failure stress in 11-direction tension              (Real>=0.0 or blank)
  double precision      :: R11Com   ! Failure stress in 11-direction compression          (Real>=0.0 or blank)
  double precision      :: R22Ten   ! Failure stress in 22-direction tension              (Real>=0.0 or blank)
  double precision      :: R22Com   ! Failure stress in 22-direction compression          (Real>=0.0 or blank)
  double precision      :: R33Ten   ! Failure stress in 33-direction tension              (Real>=0.0 or blank)
  double precision      :: R33Com   ! Failure stress in 33-direction compression          (Real>=0.0 or blank)
  double precision      :: R12Shear ! Failure stress in 12-direction shear                (Real>=0.0 or blank)
  double precision      :: R13Shear ! Failure stress in 13-direction shear                (Real>=0.0 or blank)
  double precision      :: R23Shear ! Failure stress in 23-direction shear                (Real>=0.0 or blank)              
  double precision      :: coupl12  ! Coupling coefficient 1-2                            (Real>=0.0 or blank)              
  double precision      :: coupl13  ! Coupling coefficient 1-3                            (Real>=0.0 or blank)              
  double precision      :: coupl23  ! Coupling coefficient 2-3                            (Real>=0.0 or blank)
end type fail_tsaiWu3d_type

!
! =================================================================================================
  
contains
  
  subroutine free_mem_fail(versagenskriterien)
  
    implicit none
    
    type(failurecriteria_type) :: versagenskriterien
    
    if (allocated(versagenskriterien%failTrescas)          .eq. .true.) deallocate(versagenskriterien%failTrescas)
    if (allocated(versagenskriterien%failMises)            .eq. .true.) deallocate(versagenskriterien%failMises)
    if (allocated(versagenskriterien%failMaxprincstresses) .eq. .true.) deallocate(versagenskriterien%failMaxprincstresses)
    if (allocated(versagenskriterien%failPucks)            .eq. .true.) deallocate(versagenskriterien%failPucks)
    if (allocated(versagenskriterien%failHills)            .eq. .true.) deallocate(versagenskriterien%failHills)
    if (allocated(versagenskriterien%failNorris)           .eq. .true.) deallocate(versagenskriterien%failNorris)
    if (allocated(versagenskriterien%failFibres)           .eq. .true.) deallocate(versagenskriterien%failFibres)
    if (allocated(versagenskriterien%failMaxStrains)       .eq. .true.) deallocate(versagenskriterien%failMaxStrains)
    if (allocated(versagenskriterien%failCuntzes)          .eq. .true.) deallocate(versagenskriterien%failCuntzes)
    if (allocated(versagenskriterien%failMaxStrain3Ds)     .eq. .true.) deallocate(versagenskriterien%failMaxStrain3Ds)
    if (allocated(versagenskriterien%failTsaiWu3Ds)         .eq. .true.) deallocate(versagenskriterien%failTsaiWu3Ds)
    
  end subroutine free_mem_fail

end module failurecriteria_typen
