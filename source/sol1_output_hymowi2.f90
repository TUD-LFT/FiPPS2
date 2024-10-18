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
!> Projektspezifische Ergebnisberechnung fuer HyMoWi
!
!> @details
!> Projektspezifische Ergebnisberechnung fuer die Optimierung aktiver Fluegelprofile (2D) im
!> Rahmen des HyMoWi-Projekts. Nach erfolgreicher Berechnung der Struktur unter Beruecksichtigung
!> der Fluid-Struktur-Interaktion erfolgt nachgeschaltet eine Berechnung der gewuenschten aero-
!> dynamischen Beiwerte.
!
!> @author Florian Dexl, TU Dresden, WiMi, 10.09.2018
!
!> $Id: sol1_output_hymowi2.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_hymowi2(fesim,Utot,scloop,aeroConverged)
!
! use
!
  use fesimulation_typen
  use vtk_variablen
  use failure_criteria
  
  implicit none
  
  interface
    subroutine panel2d_main(fname,readmodel,inviscid,displacements,dispfac,given_dalpha,paneling,pressures,F_total,al,ca,cd,cm,hkmax,xsep,Mcrit,output)
        character(len=*)                                                      :: fname             !< Filename of the PANEL2D input file
        logical,                                       intent(in), optional   :: readmodel         !< Flag, to read the input file
        logical,                                       intent(in), optional   :: inviscid          !< Flag, to force inviscid calculation
        double precision, dimension(:,:),              intent(in), optional   :: displacements     !< Displacements to be added on node coordinates
        double precision,                              intent(in), optional   :: dispfac           !< Factor with which displacements are added
        double precision,                              intent(in), optional   :: given_dalpha      !< Delta Alpha (overrides value from inputfile)
        logical,                                       intent(in), optional   :: paneling          !< Flag, to enable repaneling in Xfoil
        double precision, dimension(:,:), allocatable, intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        double precision, dimension(:), allocatable,   intent(out), optional  :: F_total           !< Sum of pressure times panel area
        double precision, dimension(:), allocatable,   intent(out), optional  :: al                !< Anstellwinkel
        double precision, dimension(:), allocatable,   intent(out), optional  :: ca                !< Auftriebsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cd                !< Widerstandsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cm                !< Momentenbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: hkmax             !< Maximaler kinematischer Formfaktor
        double precision, dimension(:), allocatable,   intent(out), optional  :: xsep              !< Position der turbulenten Abloesung
        double precision, dimension(:), allocatable,   intent(out), optional  :: Mcrit             !< Kritische MACH-Zahl
        logical,                                       intent(in), optional   :: output            !< Flag, if PANEL2D should write any output
    end subroutine panel2d_main

    subroutine xfoilwrapper_main(fname,readmodel,inviscid,displacements,dispfac,given_dalpha,paneling,pressures,F_total,al,ca,cd,cm,hkmax,Mcrit,output)
        character(len=*)                                                      :: fname             !< Filename of the XFOIL_WRAPPER input file
        logical,                                       intent(in), optional   :: readmodel         !< Flag, to read the input file
        logical,                                       intent(in), optional   :: inviscid          !< Flag, to force inviscid calculation
        double precision, dimension(:,:),              intent(in), optional   :: displacements     !< Displacements to be added on node coordinates
        double precision,                              intent(in), optional   :: dispfac           !< Factor with which displacements are added
        double precision,                              intent(in), optional   :: given_dalpha      !< Delta Alpha (overrides value from inputfile)
        logical,                                       intent(in), optional   :: paneling          !< Flag, to enable repaneling in Xfoil
        double precision, dimension(:,:), allocatable, intent(out), optional  :: pressures         !< Manometer pressure for all panels and cases
        double precision, dimension(:), allocatable,   intent(out), optional  :: F_total           !< Sum of pressure times panel area
        double precision, dimension(:), allocatable,   intent(out), optional  :: al                !< Anstellwinkel
        double precision, dimension(:), allocatable,   intent(out), optional  :: ca                !< Auftriebsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cd                !< Widerstandsbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: cm                !< Momentenbeiwert
        double precision, dimension(:), allocatable,   intent(out), optional  :: hkmax             !< Maximaler kinematischer Formfaktor
        double precision, dimension(:), allocatable,   intent(out), optional  :: Mcrit             !< Kritische MACH-Zahl
        logical,                                       intent(in), optional   :: output            !< Flag, if PANEL2D should write any output
    end subroutine xfoilwrapper_main
  end interface
  
  type(fe_simulation), intent(in)                         :: fesim
  double precision, dimension(fesim%num_dof), intent(in)  :: Utot
  integer,intent(in)                                      :: scloop
  logical, intent(in)                                     :: aeroConverged
  
  integer                                                 :: ii,jj
  integer                                                 :: aeroloadID
  integer                                                 :: err_code=0

  double precision                                        :: masse
  double precision, allocatable, dimension(:)             :: alpha, ca, cd, cm, hkmax, xsep, Mcrit
  double precision, dimension(3)                          :: node_deformed
  integer, parameter                                      :: outfem = 55, outdisp = 56
  integer                                                 :: step = 0
  character(len=56)                                       :: filename

  double precision                                        :: tResFac
  integer                                                 :: elem, tFailTyp
  
  double precision, dimension(:,:), allocatable           :: stress
  double precision, dimension(:), allocatable             :: nodal_temperatures
  double precision, dimension(:), allocatable             :: beam2temps
  double precision, dimension(:), allocatable             :: resFac
  double precision, dimension(:), allocatable             :: thermal_force
  double precision, dimension(:), allocatable             :: thermal_stress
  integer, dimension(:), allocatable                      :: failTyp
  double precision                                        :: RFmin
  double precision, dimension(:), allocatable             :: tstress,tstrain
  
  aeroloadID = fesim%lasten%subcases(scloop)%aeroloadID
  
  if (scloop .eq. 1) then
    OPEN(outfem, file = 'output_FEM_HyMoWi.txt', STATUS='REPLACE')

    masse = 0.d0

    do ii = 1, size(fesim%eigenschaften%pbeams,1)
        masse = masse + fesim%eigenschaften%pbeams(ii)%weight
    end do

    write(outfem,'(A21,I25)')    'max. Lastfallanzahl: ', fesim%num_subcases
    write(outfem,'(A21,E25.18)') 'Ges.-masse:          ', masse
  else
    OPEN(outfem, file = 'output_FEM_HyMoWi.txt', STATUS='OLD', POSITION='APPEND')
  end if

  if (fesim%calculateTSE .eq. .true.) then
    write(outfem,*)
    write(outfem,'(A11,E25.18)') 'TSE [-]:   ', fesim%ergebnisse%tse
  end if

  if (aeroConverged .eq. .true.) then

    if (fesim%lasten%aeroload2ds(aeroloadID)%mthd .eq. 1) then

      write(filename,'(A11,I4.4)') "panel2d_sc_", scloop
      call panel2d_main(fname=filename,readmodel=.false.,inviscid=.false.,displacements=fesim%internals%aeroDispl,paneling=.true.,al=alpha,ca=ca,cd=cd, cm=cm,hkmax=hkmax,xsep=xsep,Mcrit=Mcrit,output=fesim%ausgabe%outputVTK)
      
    else if ((fesim%lasten%aeroload2ds(aeroloadID)%mthd .eq. 2) .or. (fesim%lasten%aeroload2ds(aeroloadID)%mthd .eq. 3)) then

      write(filename,'(A11,I4.4)') "panel2d_sc_", scloop
      call xfoilwrapper_main(fname=filename,readmodel=.false.,inviscid=.false.,displacements=fesim%internals%aeroDispl,paneling=.true.,al=alpha,ca=ca,cd=cd, cm=cm,hkmax=hkmax,Mcrit=Mcrit,output=fesim%ausgabe%outputVTK)
      allocate(xsep(1))
      xsep(1) = 1.d0

    end if

    write(outfem,*)
    write(outfem,'(A11,E25.18)') 'AOA [deg]: ', alpha(1)
    write(outfem,'(A11,E25.18)') 'CA  [-]  : ', ca(1)
    write(outfem,'(A11,E25.18)') 'CD  [-]  : ', cd(1)
    write(outfem,'(A11,E25.18)') 'CM  [-]  : ', cm(1)
    write(outfem,'(A11,E25.18)') 'HK  [-]  : ', hkmax(1)
    write(outfem,'(A11,E25.18)') 'Mcrit [-]: ', Mcrit(1)
    write(outfem,'(A11,E25.18)') 'dxsep [-]: ', (1.d0 - xsep(1))
      
    deallocate(ca,cd,cm,hkmax,Mcrit)
    
    step = step + 1
    write(filename,'(A14,I4.4,A4)') "displacements_", step, ".txt"
    OPEN(outdisp, file = TRIM(filename), STATUS='UNKNOWN')
    write(outdisp,'(I10)') fesim%num_nodes
    do ii = 1,fesim%num_nodes
        node_deformed(:) = fesim%knoten%nodes(ii)%coords(:)
        node_deformed(:) = node_deformed(:) + Utot(((ii-1)*6+1):((ii-1)*6+3))
        write(outdisp,'(3E25.18)') (node_deformed(jj), jj=1,3)
    end do
    close(outdisp)
    
    allocate(nodal_temperatures(fesim%num_nodes))
    nodal_temperatures = 0.d0
  
    if (fesim%is_temperature == .true.) then
        call get_node_temperatures(fesim,fesim%lasten%subcases(scloop)%loadid, nodal_temperatures)
    end if
    
    allocate(beam2temps(size(fesim%elemente%beam2s,1)))
    beam2temps = 0.d0

    ! Get elemental temperatures at beam2 elements
    if (fesim%is_beam2temp .eqv. .true.) call get_beam2_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, beam2temps)
    
    if (fesim%is_beam2 == .true.) then
    
      allocate(stress(fesim%num_elements, 7), resFac(fesim%num_elements), failTyp(fesim%num_elements), tstress(3), tstrain(3), thermal_force(fesim%num_elements), thermal_stress(fesim%num_elements))
      
      stress = 0.d0
      resFac = 0.d0
      thermal_force  = 0.d0
      thermal_stress = 0.d0
    
      call beam2_stress(fesim, Utot(:), nodal_temperatures, beam2temps, stress, thermal_force, thermal_stress)
      
      do elem = 1, size(fesim%elemente%beam2s,1)
      
        tstrain = 0.d0
        tstress = 0.d0
        tstress(1) = stress(fesim%elemente%beam2s(elem)%eid,4)
      
        call getRF_mat1(tstress, &
                      & tstrain, &
                      & fesim%eigenschaften%pbeams(fesim%elemente%beam2s(elem)%int_pid)%intMat1ID, &
                      & resFac(fesim%elemente%beam2s(elem)%eid), &
                      & failTyp(fesim%elemente%beam2s(elem)%eid), &
                      & fesim%versagenskriterien, &
                      & fesim%materialien%mat1s)
                      
        tstress(1) = stress(fesim%elemente%beam2s(elem)%eid,5)
      
        call getRF_mat1(tstress, &
                      & tstrain, &
                      & fesim%eigenschaften%pbeams(fesim%elemente%beam2s(elem)%int_pid)%intMat1ID, &
                      & tResFac, &
                      & tFailTyp, &
                      & fesim%versagenskriterien, &
                      & fesim%materialien%mat1s)
        
        if (tResFac .LT. resFac(fesim%elemente%beam2s(elem)%eid)) then
          resFac(fesim%elemente%beam2s(elem)%eid) = tResFac
          failTyp(fesim%elemente%beam2s(elem)%eid) = tFailTyp
        end if
        
      end do
      
      write(outfem,'(A11,E25.18)') 'minRF [-]: ', minval(resFac(fesim%elemente%beam2s(1)%eid:fesim%elemente%beam2s(size(fesim%elemente%beam2s,1))%eid))
      
      deallocate(nodal_temperatures, beam2temps, stress, resFac, failTyp, tstress, tstrain, thermal_force, thermal_stress)

    end if
  end if
  
  close(outfem)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_hymowi2'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_hymowi2
