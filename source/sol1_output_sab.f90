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
!> $Id: sol1_output_sab.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine sol1_output_sab(fesim,Utot,scloop)
!
! use
!
    use globale_variablen, ONLY : shortoutFile
    use fesimulation_typen
    use vtk_variablen
    use failure_criteria
    
    implicit none
    
    type(fe_simulation), intent(inout)                     :: fesim
    double precision, dimension(fesim%num_dof), intent(in) :: Utot
    integer,intent(in)                                     :: scloop
    
    integer                                                :: ii, err_code = 0
    
    double precision                                       :: masse, tResFac
    integer                                                :: elem, shellResPos, tFailTyp
  
    double precision, dimension(:,:), allocatable          :: strain
    double precision, dimension(:,:), allocatable          :: stress
    double precision, dimension(:), allocatable            :: nodal_temperatures
    double precision, dimension(:), allocatable            :: beam2temps, quad8temps, lsolid20temps
    double precision, dimension(:), allocatable            :: resFac
    double precision, dimension(:), allocatable            :: thermal_force
    double precision, dimension(:), allocatable            :: thermal_stress
    integer, dimension(:), allocatable                     :: layNum, failTyp
    double precision, dimension(3)                         :: cent1, cent2, center
    double precision                                       :: maxXDis, maxYDis, maxZDis
    double precision                                       :: RFmin
    double precision, dimension(:), allocatable            :: tstress,tstrain
    
    if (scloop .eq. 1) then
        OPEN(shortoutFile, file = 'output_FEM.txt', STATUS='UNKNOWN')
        
        masse = 0.d0

        do ii = 1, size(fesim%eigenschaften%pcomps,1)
            masse = masse + fesim%eigenschaften%pcomps(ii)%weight
        end do

        do ii = 1, size(fesim%eigenschaften%pshells,1)
            masse = masse + fesim%eigenschaften%pshells(ii)%weight
        end do

        do ii = 1, size(fesim%eigenschaften%pbeams,1)
            masse = masse + fesim%eigenschaften%pbeams(ii)%weight
        end do

        write(shortoutFile,'(A21,E25.18)') 'Ges.-masse:          ', masse
        write(shortoutFile,'(A21,I25)')    'max. Lastfallanzahl: ', fesim%num_subcases
    else
        OPEN(shortoutFile, file = 'output_FEM.txt', STATUS='OLD', POSITION='APPEND')
    end if
    
    maxXDis = 0.d0
    maxYDis = 0.d0
    maxZDis = 0.d0
    do ii = 1, fesim%num_dof/6
      if (ABS(Utot(6*(ii-1) + 1)) > maxXDis) then
        maxXDis = ABS(Utot(6*(ii-1) + 1))
      end if
      if (ABS(Utot(6*(ii-1) + 2)) > maxYDis) then
        maxYDis = ABS(Utot(6*(ii-1) + 2))
      end if
      if (ABS(Utot(6*(ii-1) + 3)) > maxZDis) then
        maxZDis = ABS(Utot(6*(ii-1) + 3))
      end if
    end do
    
    write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', max. X-Displacement   :  ', maxXDis
    write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', max. Y-Displacement   :  ', maxYDis
    write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', max. Z-Displacement   :  ', maxZDis

    if (fesim%calculateTSE == .true.) then
      write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', total strain energy   :  ', fesim%ergebnisse%tse
    end if
    
    allocate(nodal_temperatures(fesim%num_nodes))
    nodal_temperatures = 0.d0
  
    if (fesim%is_temperature == .true.) then
        call get_node_temperatures(fesim,fesim%lasten%subcases(scloop)%loadid, nodal_temperatures)
    end if
    
    allocate(beam2temps(size(fesim%elemente%beam2s,1)), quad8temps(size(fesim%elemente%quad8s,1)), lsolid20temps(size(fesim%elemente%lsolid20s,1)))
    beam2temps = 0.d0
    quad8temps = 0.d0
    lsolid20temps = 0.d0

    ! Get elemental temperatures at beam2 elements
    if (fesim%is_beam2temp .eqv. .true.) call get_beam2_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, beam2temps)
    ! Get elemental temperatures at quad8 elements
    if (fesim%is_quad8temp .eqv. .true.) call get_quad8_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, quad8temps)
    ! Get elemental temperatures at lsolid20 elements
    if (fesim%is_lsolid20temp .eqv. .true.) call get_lsolid20_temperatures(fesim, fesim%lasten%subcases(scloop)%loadid, lsolid20temps)
    
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
      
      write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', min. B2-Reservefactor :  ', minval(resFac(fesim%elemente%beam2s(1)%eid:fesim%elemente%beam2s(size(fesim%elemente%beam2s,1))%eid))
      write(shortoutFile,'(A8,I3,A27,I25)')    'Subcase ', scloop, ', min. B2-Reservefac EID:  ', minloc(resFac(fesim%elemente%beam2s(1)%eid:fesim%elemente%beam2s(size(fesim%elemente%beam2s,1))%eid)) + fesim%elemente%beam2s(1)%eid - 1
      
      deallocate(stress, resFac, failTyp, tstress, tstrain, thermal_force, thermal_stress)

    end if

    if (fesim%is_quad8 == .true.) then
        
        ! Berechnen der Reservefaktoren für die geschichteten Elemente
        
        allocate(resFac(fesim%num_elements),layNum(fesim%num_elements),failTyp(fesim%num_elements))
        
        resFac = 0.d0
        layNum = -1
        
        call quad8_results_elem_lay(fesim, Utot(:), nodal_temperatures, quad8temps, resFac, layNum, failTyp)
        
        ! Berechnen der Reservefaktoren für die isotropen Elemente
        
        do shellResPos = 1, 2
        
            allocate(strain(fesim%num_elements, 6))
            allocate(stress(fesim%num_elements, 6))

            strain  = 0.d0
            stress  = 0.d0
        
            call quad8_results_elem(fesim, Utot(:), nodal_temperatures, quad8temps, strain, stress, shellResPos, 0, .true.)
            
            do elem = 1, size(fesim%elemente%quad8s,1)
            
                if (fesim%elemente%quad8s(elem)%propType == 1) then
                
                    call getRF_mat1(stress(fesim%elemente%quad8s(elem)%eid,1:3), &
                            & strain(fesim%elemente%quad8s(elem)%eid,1:3), &
                            & fesim%eigenschaften%pshells(fesim%elemente%quad8s(elem)%int_pid)%intMat1ID, &
                            & tResFac, &
                            & tFailTyp, &
                            & fesim%versagenskriterien, &
                            & fesim%materialien%mat1s)
                    layNum(fesim%elemente%quad8s(elem)%eid) = 0
                    
                    if (shellResPos == 1) then
                        resFac(fesim%elemente%quad8s(elem)%eid) = tResFac
                        failTyp(fesim%elemente%quad8s(elem)%eid) = tFailTyp
                    else
                        if (tResFac .LT. resFac(fesim%elemente%quad8s(elem)%eid)) then
                            resFac(fesim%elemente%quad8s(elem)%eid) = tResFac
                            failTyp(fesim%elemente%quad8s(elem)%eid) = tFailTyp
                        end if
                    end if
                        
                end if
            
            end do
            
            deallocate(stress, strain)
            
        end do
        
        do elem = 1, size(fesim%elemente%quad8s,1)
        
            cent1(1:3) = 0.5d0*(fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(1))%coords(1:3) + fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(3))%coords(1:3))
            cent2(1:3) = 0.5d0*(fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(2))%coords(1:3) + fesim%knoten%nodes(fesim%elemente%quad8s(elem)%nids(4))%coords(1:3))
            center(1:3) = 0.5d0*(cent1(1:3) + cent2(1:3))
        
            if (center(1) .LT. fesim%ausgabe%xmin .or. center(1) .GT. fesim%ausgabe%xmax) then
                resFac(fesim%elemente%quad8s(elem)%eid) = 1.d300
            end if
        
            if (center(2) .LT. fesim%ausgabe%ymin .or. center(2) .GT. fesim%ausgabe%ymax) then
                resFac(fesim%elemente%quad8s(elem)%eid) = 1.d300
            end if
        
            if (center(3) .LT. fesim%ausgabe%zmin .or. center(3) .GT. fesim%ausgabe%zmax) then
                resFac(fesim%elemente%quad8s(elem)%eid) = 1.d300
            end if
        
        end do
        
        RFmin = minval(resFac(fesim%elemente%quad8s(1)%eid:fesim%elemente%quad8s(size(fesim%elemente%quad8s,1))%eid))
        
        write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', min. Q8-Reservefactor :  ', RFmin
        write(shortoutFile,'(A8,I3,A27,I25)')    'Subcase ', scloop, ', min. Q8-Reservefac EID:  ', minloc(resFac(fesim%elemente%quad8s(1)%eid:fesim%elemente%quad8s(size(fesim%elemente%quad8s,1))%eid)) + fesim%elemente%quad8s(1)%eid - 1
        
        if (RFmin .LT. 0.9d0) then
            fesim%internals%failed = .TRUE.
        end if
        
        deallocate(resFac, layNum, failTyp)

    end if
    
    if ( fesim%is_lsolid20  == .true.) then

      allocate(resFac(fesim%num_elements),layNum(fesim%num_elements),failTyp(fesim%num_elements))
      
      resFac = 0.d0
      layNum = -1
      
      call lsolid20_results_elem_lay(fesim, Utot(:), nodal_temperatures, lsolid20temps, resFac, layNum, failTyp, scloop)
        
      write(shortoutFile,'(A8,I3,A27,E25.18)') 'Subcase ', scloop, ', min. S2-Reservefactor :  ', minval(resFac(fesim%elemente%lsolid20s(1)%eid:fesim%elemente%lsolid20s(size(fesim%elemente%lsolid20s,1))%eid))
      write(shortoutFile,'(A8,I3,A27,I25)')    'Subcase ', scloop, ', min. S2-Reservefac EID:  ', minloc(resFac(fesim%elemente%lsolid20s(1)%eid:fesim%elemente%lsolid20s(size(fesim%elemente%lsolid20s,1))%eid)) + fesim%elemente%lsolid20s(1)%eid - 1
      
      deallocate(resFac, layNum, failTyp)
      
    end if
  
    deallocate(nodal_temperatures, beam2temps, quad8temps, lsolid20temps)
  
    close(shortoutFile)
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'sol1_output_sab'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine sol1_output_sab
