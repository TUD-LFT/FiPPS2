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
SUBROUTINE xfw_write_xfoil_input(plot,reynold,mach,n_crit,xtr_up,xtr_lo,n_iter,visc,pane,autotrim,value,fname_af,fname_rs,fname_hk,fname_cp,fname,start_alpha)
  
  implicit none

  logical                , intent(in)                   :: plot, visc, pane, autotrim
  double precision       , intent(in)                   :: reynold, mach, xtr_up, xtr_lo, value
  integer                , intent(in)                   :: n_crit, n_iter
  character(len=256)     , intent(in)                   :: fname_af, fname_rs, fname_hk, fname_cp, fname  
  double precision       , intent(in),    optional      :: start_alpha

  integer                                               :: unit_xb = 21 

!--------------------------------------------------------------------------------------------------------------------------

  open(unit_xb, file = './' // TRIM(fname), status = 'replace')                     !*** Oeffnen der Xfoil-Batch-Datei

!***************************************************************************************************
!               Schreibe Xfoil-Batch-Datei
!***************************************************************************************************

  write(unit_xb,'(A4)') 'PLOP'                                                      !*** Wechseln in Menue Plotoptionen
  write(unit_xb,'(A,X,L1)') 'G', plot                                               !*** Setze Grafikoption fuer Xfoil
  write(unit_xb,*)                                                                  !*** Wechseln in OPER-Ebene
  write(unit_xb,'(A4,X,A<len_trim(fname_af)>)') 'LOAD', trim(fname_af)              !*** Laden der Profilkoordinaten
  if (pane) then
    write(unit_xb,'(A4)') 'PANE'                                                    !*** Neuvernetzen des Profils
  end if
  write(unit_xb,'(A4)') 'OPER'                                                      !*** Wechseln in Rechenmodus
  write(unit_xb,'(A2,E12.4)') 'Re', reynold                                         !*** Schreibe Reynolds-Zahl
  write(unit_xb,'(A,E12.4)') 'M', mach                                              !*** MACH-Zahl
  write(unit_xb,'(A4,X,I0)') 'ITER', n_iter                                         !*** Anzahl Iterationen
  if (visc) then
    write(unit_xb,'(A4)') 'Visc'                                                    !*** Wechseln in viskosen Modus
    write(unit_xb,'(A4)') 'VPAR'                                                    !*** Wechseln in Einstellungsmodus viskose Parameter
    write(unit_xb,'(A,X,I2)') 'N', n_crit                                           !*** Setzen der kritischen n-Zahl
    write(unit_xb,'(A3,2E12.4)') 'XTR', xtr_up, xtr_lo                              !*** Erzwinge Transition
    write(unit_xb,*)                                                                !*** Wechseln in OPER-Ebene
  end if
  if (present(start_alpha)) then                                                    !*** Wenn Start-Anstellwinkel vorgegeben
    write(unit_xb,'(A4,E12.4)') 'Alfa', start_alpha                                 !*** Gebe Anstellwinkel vor
  end if
  write(unit_xb,'(A4)') 'Pacc'                                                      !*** Automatische Akkumulierung der aktiven Polare
  write(unit_xb,'(A<len_trim(fname_rs)>)') trim(fname_rs)                           !*** Definition der Ausgabedatei
  write(unit_xb,*)                                                                  !*** Keine Dump-Datei
  if (autotrim .EQV. .TRUE.) then
    write(unit_xb,'(A2,E12.4)') 'Cl', value                                         !*** Gebe Auftriebsbeiwert vor
  else
    write(unit_xb,'(A4,E12.4)') 'Alfa', value                                       !*** Gebe Anstellwinkel vor
  end if
  write(unit_xb,'(A4,X,A<len_trim(fname_cp)>)') 'CPWR', trim(fname_cp)              !*** Schreibe Druckverteilung cp(x) in Datei
  if (visc) then
    write(unit_xb,'(A4)') 'VPlo'                                                      !*** Wechseln in Menue fuer Grenzschichtgroessen
    write(unit_xb,'(A)') 'H'                                                          !*** Plotte Formfaktor H
    write(unit_xb,'(A4,X,A<len_trim(fname_hk)>)') 'DUMP', trim(fname_hk)              !*** Schreibe Formfaktor H in Datei
    write(unit_xb,*)                                                                  !*** Wechseln in OPER-Ebene
  end if
  write(unit_xb,*)                                                                  !*** Wechseln in XFOIL-Ebene
  write(unit_xb,'(A4)') 'QUIT'                                                      !*** Xfoil beenden

!***************************************************************************************************
!               Schliesse Datei
!***************************************************************************************************

  close(unit_xb)                                                                    !*** Schliessen der Datei

END SUBROUTINE xfw_write_xfoil_input
