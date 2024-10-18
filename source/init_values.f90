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
!> Subroutine initializes neccesary values for the calculation
!
!> @details
!> The subroutine computes and initializes some global variables of the
!> actual problem
!> Revision:    Martin Rädel   21.07.2010
!>              TU Dresden, WiMi
!>              -> Überschreiben fesim%eigenschaften%pcomps(ii)%lay aus lam8-Karten eingebaut
!
!> @author Martin Rädel, TU Dresden, wissenschaftlicher Mitarbeiter 21.07.2010
!
!> $Id: init_values.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine init_values(fesim)
!
! Use
!
  use konstanten
  use fesimulation_typen
!
! =================================================================================================
!
  implicit none
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
! Input+Output
!
  type(fe_simulation)  :: fesim
!
! Internal
!
  integer :: ii
  integer :: oldCPSet
  integer :: lam8id, nrlay, jj
  integer :: lam20id
  logical :: flag
  
  integer       :: err_code=0
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
  ! Compute number of dof's in whole model
  !
  ! Remark: all nodes are supposed to have 6 degrees of freedom
  !
  fesim%num_dof     = ndof*size(fesim%knoten%nodes,1)

  fesim%num_elements = 0
  if (allocated(fesim%elemente%beam2s) .eq. .true.) fesim%num_elements = fesim%num_elements + size(fesim%elemente%beam2s,1)
  if (allocated(fesim%elemente%quad8s) .eq. .true.) fesim%num_elements = fesim%num_elements + size(fesim%elemente%quad8s,1)
  if (allocated(fesim%elemente%lsolid20s) .eq. .true.) fesim%num_elements = fesim%num_elements + size(fesim%elemente%lsolid20s,1)
  
  !
  ! Determine the number of subcases
  !
  if (fesim%is_subcase == .true.) then
    fesim%num_subcases=1
    if (size(fesim%lasten%subcases,1) > 1) then
      do ii=2,size(fesim%lasten%subcases,1)
        if (fesim%lasten%subcases(ii)%scid /= fesim%lasten%subcases(ii-1)%scid) then
          fesim%num_subcases=fesim%num_subcases+1
        end if
      end do
    end if
  else
      write(*,*) 'No subcases definded, exit program'
      err_code=1
      goto 9999
  end if
  
  fesim%lasten%subcases(1)%upmats = .true.
  do ii = 2, fesim%num_subcases
    fesim%lasten%subcases(ii)%upmats = .true.
    if (fesim%lasten%subcases(ii)%upgeom .eq. .false. .AND. &
        fesim%lasten%subcases(ii-1)%spcaddid .eq. fesim%lasten%subcases(ii)%spcaddid .AND. &
        fesim%lasten%subcases(ii-1)%mpcaddid .eq. fesim%lasten%subcases(ii)%mpcaddid) then
        fesim%lasten%subcases(ii)%upmats = .false.
    end if
  end do
  
  ! 
  ! Wenn die Reaktionskräfte bestimmt werden sollen, wird nach jedem Subcase die
  ! vollständige Steifikeitsmatrix (ohne Verschiebungsrandbedingungen) aufgebaut. 
  ! Damit muss bei jedem neuen Subcase die Matrix mit Randbedingungen neu 
  ! berechnet werden.
  !
  if (fesim%calculateReactForce .eq. .true.) then
    do ii = 2, fesim%num_subcases
      fesim%lasten%subcases(ii)%upmats = .true.
    end do
  end if
  
  !
  ! Determine the number of loadcases
  !
  fesim%num_loadcases=1                                             ! 1st row
  if (size(fesim%lasten%loads,1) > 1) then
    do ii=2,size(fesim%lasten%loads,1)
        if (fesim%lasten%loads(ii)%lcid /= fesim%lasten%loads(ii-1)%lcid) then
            flag=.false.
            do jj=1,(ii-1)
                if (fesim%lasten%loads(ii)%lcid == fesim%lasten%loads(jj)%lcid) then
                    flag=.true.
                end if
            end do
            if (flag == .false.) then
                fesim%num_loadcases = fesim%num_loadcases+1
            end if
        end if
    end do
  end if
  
  !
  ! Determine highest loadID
  !
  if (fesim%is_force .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%forces(1:size(fesim%lasten%forces,1))%lid))
  end if
  if (fesim%is_moment .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%moments(1:size(fesim%lasten%moments,1))%lid))
  end if
  if (fesim%is_p2load .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%p2loads(1:size(fesim%lasten%p2loads,1))%lid))
  end if
  if (fesim%is_p8load .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%p8loads(1:size(fesim%lasten%p8loads,1))%lid))
  end if
  if (fesim%is_p20load .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%p20loads(1:size(fesim%lasten%p20loads,1))%lid))
  end if
  if (fesim%is_aeroload2d .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%aeroload2ds(1:size(fesim%lasten%aeroload2ds,1))%lid))
  end if
  if (fesim%is_aeroload3d .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%aeroload3ds(1:size(fesim%lasten%aeroload3ds,1))%lid))
  end if
  if (fesim%is_temperature .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%nodetemps(1:size(fesim%lasten%nodetemps,1))%lid))
  end if
  if (fesim%is_beam2temp .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%beam2temps(1:size(fesim%lasten%beam2temps,1))%lid))
  end if
  if (fesim%is_quad8temp .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%quad8temps(1:size(fesim%lasten%quad8temps,1))%lid))
  end if
  if (fesim%is_lsolid20temp .eq. .true.) then
    fesim%internals%highestLID = max(fesim%internals%highestLID, maxval(fesim%lasten%lsolid20temps(1:size(fesim%lasten%lsolid20temps,1))%lid))
  end if
  
  !
  ! Determine the number of coupling sets
  !
  if (fesim%is_coupling .EQ. .true.) then
    oldCPSet = fesim%randbedingungen%couplings(1)%cpsid
    fesim%num_cps = 1
    do ii = 2, size(fesim%randbedingungen%couplings,1)
      if (fesim%randbedingungen%couplings(ii)%cpsid .ne. oldCPSet) then
        fesim%num_cps = fesim%num_cps + 1
      end if  
    end do
  end if
  
  !
  ! Calculate number of layers in pcomp from lam8 and overwrite value in pcomp
  !
  
  if (fesim%is_lam8 == .true. .AND. fesim%is_pcomp == .true.) then
    do ii=1,size(fesim%eigenschaften%pcomps)
      ! get lam8-number
      lam8id = fesim%eigenschaften%pcomps(ii)%lamid
      ! get number of layers in lam8
      nrlay=0
      do jj=1,size(fesim%laminate%lam8s)
        if (fesim%laminate%lam8s(jj)%lamid == lam8id) then
          nrlay = nrlay + 1
        end if
      end do
      ! overwrite value in pcomp card
      fesim%eigenschaften%pcomps(ii)%lay = nrlay
    end do
  end if
  
  !
  ! Calculate number of layers in fesim%eigenschaften%plsolids from lam20 and overwrite value in fesim%eigenschaften%plsolids
  !

  if (fesim%is_lam20 == .true. .AND. fesim%is_plsolid == .true.) then
    do ii=1,size(fesim%eigenschaften%plsolids)
      ! get lam20-number
      lam20id = fesim%eigenschaften%plsolids(ii)%lamid
      ! get number of layers in lam20
      nrlay=0
      do jj=1,size(fesim%laminate%lam20s)
        if (fesim%laminate%lam20s(jj)%lamid == lam20id) then
          nrlay = nrlay + 1
        end if
      end do
      ! overwrite value in plsolid card
      fesim%eigenschaften%plsolids(ii)%lay = nrlay
    end do     
  end if
  
  ! Calculate Properties for Quad8-element (ABD-matrix, surface)
  
  ! Init Properties
  if (fesim%is_quad8 == .true.) then
  
    call quad8_prepareProperties(fesim)
    
  end if
  
  if (fesim%is_beam2 == .true.) then
  
    call beam2_prepareProperties(fesim)
  
  end if
  
  if (fesim%is_lsolid20 == .true.) then
    
    call lsolid20_prepareProperties(fesim,1)

  end if
  
  if (fesim%is_failure == .true.) then
  
    call failure_prepare(fesim)
  
  end if

  if ((fesim%is_aeroload2d == .true.) .or. (fesim%is_aeroload3d == .true.)) then

    call aeroload_prepare(fesim)

  end if

!
! =================================================================================================
!
! Error handling
!  
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'get_Utot'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if
  
end subroutine init_values
