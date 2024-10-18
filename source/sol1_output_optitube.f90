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
!> $Id: sol1_output_optitube.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_optitube(Utot)

  use globale_variablen
  use vtk_variablen
  use netz_variablen
  
  implicit none
  
  integer                                            :: err_code=0
  
  integer                                            :: ii, jj
  
  double precision                                   :: y_min, y_max, eps_krit
  double precision                                   :: epstemp, epsmax, seqvtemp, seqv_max, angle
  double precision, dimension(3)                     :: epsTrans
  double precision, dimension(3,3)                   :: transMat
  double precision, dimension(num_dof, num_subcases) :: Utot
  double precision, dimension(num_elements, 6)       :: epsilon
  double precision, dimension(num_elements, 6)       :: sigma
  double precision, dimension(num_elements, 4)       :: epsilonRes
  integer, parameter                                 :: outfem = 55
  integer                                            :: mm, elem
  double precision                                   :: masse
  double precision                                   :: RFmin
  double precision, dimension(3)                     :: cent1, cent2, center
  double precision, dimension(:), allocatable        :: nodal_temperatures
  double precision, dimension(:), allocatable        :: elemental_temperatures
  
  ! Lese Eingabedatei mit y-Bereich fuer Dehnungsauswertung und kritischer Dehnung ein
  CALL input_optitube(y_min,y_max,eps_krit)
  
  allocate(nodal_temperatures(num_nodes))
  nodal_temperatures = 0.d0
  
  allocate(elemental_temperatures(num_elements))
  elemental_temperatures = 0.d0

  OPEN(outfem, file = 'output_optitube_FEM.txt', STATUS='UNKNOWN')

  masse = 0.d0

  do ii = 1, size(pcomps,1)
    masse = masse + pcomps(ii)%weight
  end do

  do ii = 1, size(pshells,1)
    masse = masse + pshells(ii)%weight
  end do

  write(outfem,'(A34,E25.18)') '#GESMAS Ges.-masse              : ', masse

  angle = 45.0d0/180.0d0*pi
  
  do ii = 1, num_subcases
  
     if (is_temperature == .true.) then
       call get_node_temperatures(subcases(ii)%loadid, nodal_temperatures)
     end if
  
     if (is_elementtemp == .true.) then
       call get_element_temperatures(subcases(ii)%loadid, elemental_temperatures)
     end if
 
    if (is_quad8 == .true.) then
    
     epsmax = 0.d0
     seqv_max = 0.d0
     
     do mm = 1,2
       epsilon = 0.d0
       sigma = 0.d0
   
       call quad8_results_elem(Utot(:,ii), nodal_temperatures, elemental_temperatures, epsilon, sigma, mm, 0, .false.)

       do elem = 1, size(quad8s,1)
       
         cent1(1:3) = 0.5d0*(nodes(quad8s(elem)%nids(1))%coords(1:3) + nodes(quad8s(elem)%nids(3))%coords(1:3))
         cent2(1:3) = 0.5d0*(nodes(quad8s(elem)%nids(2))%coords(1:3) + nodes(quad8s(elem)%nids(4))%coords(1:3))
         center(1:3) = 0.5d0*(cent1(1:3) + cent2(1:3))
       
         if (center(2) .LT. y_min .or. center(2) .GT. y_max) then
           epsilon(quad8s(elem)%eid,1:6) = 0.d0
           sigma(quad8s(elem)%eid,1:6) = 0.d0
         end if
         
         ! Berechne Dehnungen unter 45 grad         
         transMat(1,1) = dcos(angle)**2.0d0
         transMat(1,2) = dsin(angle)**2.0d0
         transMat(1,3) = dsin(angle)*cos(angle)
         transMat(2,1) = dsin(angle)**2.0d0
         transMat(2,2) = dcos(angle)**2.0d0
         transMat(2,3) = -1.0d0*dsin(angle)*dcos(angle)
         transMat(3,1) = -2.0d0*dsin(angle)*dcos(angle)
         transMat(3,2) = 2.0d0*dsin(angle)*dcos(angle)
         transMat(3,3) = dcos(angle)**2.0d0-dsin(angle)**2.0d0
         
         epsTrans(:) = matmul(transMat,epsilon(quad8s(elem)%eid,1:3))
         
         ! Dehnung unter 0 und 90 grad
         epsilonRes(quad8s(elem)%eid,1:2) = epsilon(quad8s(elem)%eid,1:2)
         ! Dehnung unter +45 und -45 grad
         epsilonRes(quad8s(elem)%eid,3:4) = epsTrans(1:2)

       end do
       
       ! Bestimme maximale Dehnung
       do jj = 1,4
         epstemp = maxval(abs(epsilonRes(:,jj)))
         if (epstemp .GT. epsmax) epsmax = epstemp
       end do
       
       seqvtemp = maxval(abs(sigma(:,6)))
       
       if (seqvtemp .GT. seqv_max) seqv_max = seqvtemp

     end do

     write(outfem,'(A10,I2.2)') '#LOADCASE_', ii
     
     write(outfem,'(A34,E25.18)') '#MAXEPS Dehn.-Max. (Elem.werte) : ', epsmax
     
     write(outfem,'(A34,E25.18)') '#MAXSEQ Span.-Max. (Elem.werte) : ', seqv_max
     
     RFmin = eps_krit/epsmax
     
     write(outfem,'(A34,E25.18)') '#MINRF min. RF                  : ', RFmin
     
     if (RFmin .LT. 0.9d0) failed = .TRUE.

    end if
    
  end do
  
  close(outfem)
  
  deallocate(nodal_temperatures, elemental_temperatures)

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_optitube'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

contains
  subroutine input_optitube(y_min,y_max,eps_krit)
  !
  ! =================================================================================================
  !
  implicit none
  !
  ! =================================================================================================
  !
  intent(out)      :: y_min, y_max, eps_krit
  double precision :: y_min, y_max, eps_krit
  character(50)    :: fname                                  ! filename with maximum length of (XX) characters
  integer          :: infem = 60
  integer          :: io_error, read_error                   ! Error-handling parameters
  integer          :: rowcount
  integer          :: ii
  integer          :: err_code=0
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
  
  fname='optitube.fipps'
  open (unit=infem,file=fname,status='old',action='read', iostat=io_error)
  if (io_error == 0) then
      read (infem,*,iostat=read_error) rowcount
      if (read_error /= 0) then
          write(*,*) 'Error reading file ', fname
          write(*,*) 'There was a problem reading the number of rows of the file'
          err_code=1
          goto 9999
      end if
      ! read values
      do ii=1,rowcount
          read (infem,'(3E23.16)',iostat=read_error) y_min, y_max, eps_krit
          if (read_error /= 0) then
              write(*,*) 'Error reading file ', fname
              write(*,*) 'There was an error processing line', ii
              err_code=2
              goto 9999
          end if
      end do
      close (infem, iostat=io_error)
      if (io_error /= 0) then
          write(*,*) 'Error closing file ', fname
          err_code=1
          goto 9999
  end if
  else
      write (*,*) 'Error opening file ', fname
      err_code=1
      goto 9999
  end if
  
  !
  ! =================================================================================================
  !
  ! Error handling
  !
  
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)                      'An error occured in subroutine'
     write(*,*)                      'input_optitube'
     write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
     write(*,*)                      'exit program '
     stop
     
  end if
  
  return
  
  end subroutine input_optitube

end subroutine sol1_output_optitube
