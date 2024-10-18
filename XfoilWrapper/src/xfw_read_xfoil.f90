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
subroutine xfw_read_xfoil(alpha,CL,CD,CM,max_HK,cpx,mach_crit,success,pane,visc,fname_rs,fname_hk,fname_cp)

  use xfw_netz_variablen
  use xfw_functions
  
  implicit none
  
  intent(in)                                                           :: pane, visc, fname_rs, fname_hk, fname_cp
  intent(out)                                                          :: alpha, CL, CD, CM, max_HK, cpx, mach_crit, success
  integer                                                              :: ii
  integer                                                              :: rstat, side
  double precision                                                     :: alpha, CL, CD, CM, max_HK, mach_crit, dump, temp
  double precision       ,dimension(n_nodes)                           :: cpx
  double precision       ,dimension(:), allocatable                    :: cp_all
  character(len=256)                                                   :: fname_rs, fname_hk, fname_cp, line
  logical                                                              :: cont_reading, success, pane, visc
  integer                                                              :: n_cp
  integer                                                              :: io_error
  integer                                                              :: err_code = 0
  integer                                                              :: unit_res = 23

!***************************************************************************************************
!               Lese Beiwerte
!***************************************************************************************************

  alpha    = huge(1.d0)
  CL       = huge(1.d0)
  CD       = huge(1.d0)
  CM       = huge(1.d0)
  max_HK   = huge(1.d0)

  ! Oeffne Datei mit Ergebnissen
  open(unit=unit_res,file=trim(fname_rs),status='old',action='read',iostat=io_error)

  ! Lese Auftriebs-, Widerstands- und Momentenbeiwert ein
  if (io_error .EQ. 0) then
    cont_reading = .true.
    success = .true.
    
    do while (cont_reading)
      read(unit_res,*) line
      line = trim(line)
      if (line(1:5) == 'alpha') cont_reading = .false.
    end do

    read(unit_res,*)

    if (visc) then
      read(unit_res,'(X,3E26.18,13X,E26.18)',iostat=rstat) alpha, CL, CD, CM
      if (rstat .NE. 0) then
        success = .false.
      end if
    else
      read(unit_res,'(X,2E26.18,26X,E13.5,E26.18)',iostat=rstat) alpha, CL, CD, CM
      CD = -1.0d0*CD
      if (rstat .NE. 0) then
        success = .false.
      end if
    end if

    ! Schliesse Datei
    close(unit_res)
  else
    success = .false.
  end if

  ! Lese maximalen Formfaktor HK der Grenzschichttheorie ein
  if ((visc .eqv. .true.) .and. (success .EQV. .true.)) then
    ! Oeffne Datei mit Ergebnissen
    open(unit=unit_res,file=trim(fname_hk),status='old',action='read',iostat=io_error)

    if (io_error .EQ. 0) then
      cont_reading = .true.
      success = .true.
      do while (cont_reading)
        read(unit_res,'(A)') line
        line = TRIM(line)
        if (line(1:6) == '#    x') cont_reading = .false.
      end do
        
      max_HK = 0.d0
      side = 1
      do while (side .LE. 2)
        read(unit_res,'(A)') line
        line = trim(line)
        if (len_trim(line) == 0) then
          side = side + 1
          cycle
        end if
        read(line,'(X,2G14.6)',iostat=rstat) dump, temp
        if (rstat .ne. 0) then
          success = .false.
        end if
        if (temp .gt. max_HK) max_HK = temp
      end do

      ! Schliesse Datei
      close(unit_res)
    else
      success = .false.
    end if
  end if

  ! Lese Druckbeiwerte ueber x ein
  if (pane .eqv. .false.) then
    ! Keine Neuvernetzung => Anzahl Druckbeiwerte muss Anzahl der Panels entsprechen
    if (success .eqv. .true.) then                             
      ! Oeffne Datei mit Ergebnissen
      open(unit=unit_res,file=trim(fname_cp),status='old',action='read',iostat=io_error)
    
      if (io_error .eq. 0) then
        read(unit_res,*)
          
        do ii = 1, n_nodes
          read(unit_res,'(A)') line
          line = trim(line)
          read(line,'(X,2F11.5)',iostat=rstat) dump, cpx(ii)
          if (rstat .NE. 0) then
            success = .false.
          end if
        end do

        n_cp = n_nodes
        allocate(cp_all(n_cp))
        cp_all(:) = cpx(:)

        ! Schliesse Datei
        close(unit_res)
      else
        success = .false.
      end if
    end if
  else
    ! Neuvernetzung => Anzahl der Panels zunaechst unbekannt
    ! Ermittle zuerst Panelanzahl und lese dann cp-Werte ein
    if (success .eqv. .true.) then                             
      ! Oeffne Datei mit Ergebnissen
      open(unit=unit_res,file=trim(fname_cp),status='old',action='read',iostat=io_error)
    
      if (io_error .eq. 0) then
        read(unit_res,*)

        n_cp = 0
        do while (.true.)
          read(unit_res,'(A)',iostat=rstat) line
          if (rstat .NE. 0) then
            exit
          else
            n_cp = n_cp + 1
          end if
        end do

        rewind(unit_res)
        read(unit_res,*)
        allocate(cp_all(n_cp))

        do ii = 1, n_cp
          read(unit_res,'(A)') line
          line = trim(line)
          read(line,'(X,2F11.5)',iostat=rstat) dump, cp_all(ii)
          if (rstat .NE. 0) then
            success = .false.
          end if
        end do
    
        ! Schliesse Datei
        close(unit_res)
      else
        success = .false.
      end if
    end if
  end if

  if (success .eqv. .true.) then
    call xfw_ma_crit(n_cp, cp_all(:), Ma, mach_crit)
  end if
  
  if (allocated(cp_all)) deallocate(cp_all)
  
  !***************************************************************************************************
  !               FEHLERBEHANDLUNG
  !***************************************************************************************************
  
  9999 continue
  
  if (err_code /= 0) then
     
     write(*,*)              'Ein Fehler trat in folgender Subroutine auf:'
     write(*,*)              'xfw_read_xfoil.'
     write(*,*)              'beende Programm'
     stop
     
  end if

end subroutine xfw_read_xfoil
