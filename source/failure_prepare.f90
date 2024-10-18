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
!> Preperation of failure criteria 
!
!> @details
!> Preperation of the failure criteria. Determine the internal failure IDs and failure type
!
!> @author Andreas Hauffe, TU Dresden, wiss. Mitarbeiter, 25.09.2017 
!
!> $Id: failure_prepare.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine failure_prepare(fesim)
! =================================================================================================
!
! Use
!
  use fesimulation_typen
  
  type(fe_simulation) :: fesim
  
  logical :: found
  
  integer :: ii,jj,mm
  integer :: err_code 
  
  err_code = 0
  
!--------------------------------------------------------------------------------------------------
!
! Mat1
!
  do ii = 1, size(fesim%materialien%mat1s,1)
  
    fesim%materialien%mat1s(ii)%ifid = -1
    fesim%materialien%mat1s(ii)%iftype = -1
    
    do mm = 1,4
    
      found = .false.
      
      if (fesim%materialien%mat1s(ii)%fid(mm) .eq. 0) then
        exit
      end if
      
      ! check Tesca
      
      if (allocated(fesim%versagenskriterien%failTrescas) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failTrescas,1)
          if (fesim%versagenskriterien%failTrescas(jj)%fid .eq. fesim%materialien%mat1s(ii)%fid(mm)) then
            fesim%materialien%mat1s(ii)%ifid(mm) = jj
            fesim%materialien%mat1s(ii)%iftype(mm) = 1
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von Mises
      
      if (allocated(fesim%versagenskriterien%failMises) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failMises,1)
          if (fesim%versagenskriterien%failMises(jj)%fid .eq. fesim%materialien%mat1s(ii)%fid(mm)) then
            fesim%materialien%mat1s(ii)%ifid(mm) = jj
            fesim%materialien%mat1s(ii)%iftype(mm) = 2
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von maximum principal Stresses
      
      if (allocated(fesim%versagenskriterien%failMaxprincstresses) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failMaxprincstresses,1)
          if (fesim%versagenskriterien%failMaxprincstresses(jj)%fid .eq. fesim%materialien%mat1s(ii)%fid(mm)) then
            fesim%materialien%mat1s(ii)%ifid(mm) = jj
            fesim%materialien%mat1s(ii)%iftype(mm) = 3
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .false.) then
        write(*,*) 'No failure criterion found for mat1 mid:', fesim%materialien%mat1s(ii)%mid
        err_code = 1
        goto 9999
      end if
  
    end do
    
  end do
  
!--------------------------------------------------------------------------------------------------
!
! Mat8
!
  do ii = 1, size(fesim%materialien%mat8s,1)
  
    fesim%materialien%mat8s(ii)%ifid = -1
    fesim%materialien%mat8s(ii)%iftype = -1
  
    do mm = 1,4
    
      found = .false.
      
      if (fesim%materialien%mat8s(ii)%fid(mm) .eq. 0) then
        exit
      end if
      
      ! check Puck
      
      if (allocated(fesim%versagenskriterien%failPucks) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failPucks,1)
          if (fesim%versagenskriterien%failPucks(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 1
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check Tesca
      
      if (allocated(fesim%versagenskriterien%failTrescas) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failTrescas,1)
          if (fesim%versagenskriterien%failTrescas(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 2
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von Mises
      
      if (allocated(fesim%versagenskriterien%failMises) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failMises,1)
          if (fesim%versagenskriterien%failMises(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 3
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von maximum principal Stresses
      
      if (allocated(fesim%versagenskriterien%failMaxprincstresses) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failMaxprincstresses,1)
          if (fesim%versagenskriterien%failMaxprincstresses(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 4
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von Hill
      
      if (allocated(fesim%versagenskriterien%failHills) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failHills,1)
          if (fesim%versagenskriterien%failHills(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 5
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von Norris
      
      if (allocated(fesim%versagenskriterien%failNorris) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failNorris,1)
          if (fesim%versagenskriterien%failNorris(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 6
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von Fibre failure
      
      if (allocated(fesim%versagenskriterien%failFibres) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failFibres,1)
          if (fesim%versagenskriterien%failFibres(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 7
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von Maximum strain
      
      if (allocated(fesim%versagenskriterien%failMaxStrains) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failMaxStrains,1)
          if (fesim%versagenskriterien%failMaxStrains(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 8
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check von 2D Cuntze
      
      if (allocated(fesim%versagenskriterien%failCuntzes) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failCuntzes,1)
          if (fesim%versagenskriterien%failCuntzes(jj)%fid .eq. fesim%materialien%mat8s(ii)%fid(mm)) then
            fesim%materialien%mat8s(ii)%ifid(mm) = jj
            fesim%materialien%mat8s(ii)%iftype(mm) = 9
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .false.) then
        write(*,*) 'No failure criterion found for mat8 mid:', fesim%materialien%mat8s(ii)%mid
        err_code = 1
        goto 9999
      end if
  
    end do
    
  end do

!--------------------------------------------------------------------------------------------------
!
! Mat20
!
  do ii = 1, size(fesim%materialien%mat20s,1)
  
    fesim%materialien%mat20s(ii)%ifid = -1
    fesim%materialien%mat20s(ii)%iftype = -1
  
    do mm = 1,4
    
      found = .false.
      
      if (fesim%materialien%mat20s(ii)%fid(mm) .eq. 0) then
        exit
      end if
      
      ! check three dimensional maximum strain criterion
      
      if (allocated(fesim%versagenskriterien%failMaxStrain3Ds) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failMaxStrain3Ds,1)
          if (fesim%versagenskriterien%failMaxStrain3Ds(jj)%fid .eq. fesim%materialien%mat20s(ii)%fid(mm)) then
            fesim%materialien%mat20s(ii)%ifid(mm) = jj
            fesim%materialien%mat20s(ii)%iftype(mm) = 108
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      ! check three dimensional Tsai-Wu-criterion
      
      if (allocated(fesim%versagenskriterien%failTsaiWu3Ds) .eq. .true.) then
        do jj = 1,size(fesim%versagenskriterien%failTsaiWu3Ds,1)
          if (fesim%versagenskriterien%failTsaiWu3Ds(jj)%fid .eq. fesim%materialien%mat20s(ii)%fid(mm)) then
            fesim%materialien%mat20s(ii)%ifid(mm) = jj
            fesim%materialien%mat20s(ii)%iftype(mm) = 110
            found = .true.
            exit
          end if
        end do
      end if
      
      if (found .eq. .true.) then
        cycle
      end if
      
      if (found .eq. .false.) then
        write(*,*) 'No failure criterion found for mat20 mid:', fesim%materialien%mat20s(ii)%mid
        err_code = 1
        goto 9999
      end if
  
    end do
    
  end do
!
! =================================================================================================
!
! Error handling
!
9999 continue

  if (err_code /= 0) then
   
    write(*,*)                      'An error occured in subroutine'
    write(*,*)                      'failure_prepare'
    write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
    write(*,*)                      'exit program '
    stop
   
  end if

end subroutine failure_prepare
