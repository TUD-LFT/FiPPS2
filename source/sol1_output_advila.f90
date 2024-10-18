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
!> Projektspezifische Ergebnisberechnung fuer AdviLa
!
!> @details
!> Projektspezifische Ergebnisberechnung fuer AdviLa
!
!> @author Florian Dexl, TU Dresden, WiMi, 02.07.2020
!
!> $Id: sol1_output_hymowi2.f90 423 2020-03-26 16:43:25Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 423 $
!> $Date: 2020-03-26 17:43:25 +0100 (Do, 26. Mär 2020) $
!
! =================================================================================================
subroutine sol1_output_advila(fesim,Utot,scloop)
!
! use
!
  use fesimulation_typen
  use vtk_variablen
  use failure_criteria
  
  implicit none
  
  type(fe_simulation), intent(in)                         :: fesim
  double precision, dimension(fesim%num_dof), intent(in)  :: Utot
  integer,intent(in)                                      :: scloop
  
  integer                                                 :: ii
  integer                                                 :: rowcount
  integer                                                 :: io_error, read_error
  integer                                                 :: err_code=0

  integer, dimension(:), allocatable                      :: target_elements
  double precision, dimension(:), allocatable             :: quad8temps, nodal_temperatures, tseElementalVec

  double precision                                        :: tse_sum
  integer, parameter                                      :: infem= 56, outfem = 55
  character(len=56)                                       :: filename

  allocate(quad8temps(size(fesim%elemente%quad8s,1)))
  quad8temps = 0.d0
  ! Get elemental temperatures at quad8 elements
  if (fesim%is_quad8temp .eqv. .true.) call get_quad8_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, quad8temps)
  
  allocate(nodal_temperatures(fesim%num_nodes))
  nodal_temperatures = 0.d0
  if (fesim%is_temperature == .true.) call get_node_temperatures(fesim,fesim%lasten%subcases(scloop)%loadid, nodal_temperatures) 

  allocate(tseElementalVec(fesim%num_elements))
  tseElementalVec = 0.d0
  call quad8_strainEnergy_control(fesim, Utot(:), nodal_temperatures, quad8temps, tseElementalVec)

  deallocate(quad8temps, nodal_temperatures)

  if (scloop .eq. 1) then
    OPEN(outfem, file = 'output_FEM_AdviLa.txt', STATUS='REPLACE')
  else
    OPEN(outfem, file = 'output_FEM_AdviLa.txt', STATUS='OLD', POSITION='APPEND')
  end if

  filename='advila_tseElements_all.fipps'
  open (unit=infem,file=filename,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
      read (infem,*,iostat=read_error) rowcount
      if (read_error /= 0) then
          write(*,*) 'Error reading file ', filename
          write(*,*) 'There was a problem reading the number of rows of the file'
          err_code=1
          goto 9999
      end if
      ! read values
      allocate(target_elements(rowcount))
      do ii=1,rowcount
          read (infem,'(I10)',iostat=read_error) target_elements(ii)
          if (read_error /= 0) then
              write(*,*) 'Error reading file ', filename
              write(*,*) 'There was an error processing line', ii
              err_code=2
              goto 9999
          end if
      end do
      close (infem, iostat=io_error)
      if (io_error /= 0) then
          write(*,*) 'Error closing file ', filename
          err_code=1
          goto 9999
  end if
  else
      write (*,*) 'Error opening file ', filename
      err_code=1
      goto 9999
  end if
  
  tse_sum = 0.d0
  do ii = 1, size(target_elements)
    tse_sum = tse_sum + tseElementalVec(target_elements(ii))
  end do
  write(outfem,'(A8,I3,A15,E25.18)')    'Subcase ', scloop, ', TSE - all:   ', tse_sum

  deallocate(target_elements)
  
  filename='advila_tseElements_inner.fipps'
  open (unit=infem,file=filename,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
      read (infem,*,iostat=read_error) rowcount
      if (read_error /= 0) then
          write(*,*) 'Error reading file ', filename
          write(*,*) 'There was a problem reading the number of rows of the file'
          err_code=1
          goto 9999
      end if
      ! read values
      allocate(target_elements(rowcount))
      do ii=1,rowcount
          read (infem,'(I10)',iostat=read_error) target_elements(ii)
          if (read_error /= 0) then
              write(*,*) 'Error reading file ', filename
              write(*,*) 'There was an error processing line', ii
              err_code=2
              goto 9999
          end if
      end do
      close (infem, iostat=io_error)
      if (io_error /= 0) then
          write(*,*) 'Error closing file ', filename
          err_code=1
          goto 9999
  end if
  else
      write (*,*) 'Error opening file ', filename
      err_code=1
      goto 9999
  end if

  tse_sum = 0.d0
  do ii = 1, size(target_elements)
    tse_sum = tse_sum + tseElementalVec(target_elements(ii))
  end do
  write(outfem,'(A8,I3,A15,E25.18)')    'Subcase ', scloop, ', TSE - inner: ', tse_sum

  deallocate(target_elements)
  
  deallocate(tseElementalVec)
  
  close(outfem)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_advila'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_advila
