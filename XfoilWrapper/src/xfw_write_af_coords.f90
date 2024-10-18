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
SUBROUTINE xfw_write_af_coords(fname_af)

  use xfw_netz_variablen

  implicit none
  
  character(len=256), intent(in)                        :: fname_af
  integer                                               :: ii, jj
  integer                                               :: unit_af = 22

!--------------------------------------------------------------------------------------------------------------------------

!***************************************************************************************************
!               Schreibe Profilkoordinaten
!***************************************************************************************************

  ! Oeffnen der Profilkoordinaten-Datei
  open(unit_af,file=trim(fname_af),status = 'replace')

  ! Schreibe Profilname
  write(unit_af,'(A12)') 'XfoilWrapper'
  
  do ii=1,n_nodes
      write(unit_af,'(2E24.16)') (nodes(ii,jj), jj=1,2)
  end do 

  ! Schliessen der Datei
  close(unit_af)                                                                    

END SUBROUTINE xfw_write_af_coords
