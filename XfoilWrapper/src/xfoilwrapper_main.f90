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
subroutine xfoilwrapper_main(input_file, readmodel, inviscid, displacements, dispfac, given_dalpha, paneling, pressures, F_total, al, ca, cd, cm, hkmax, Macrit, output)

! Module laden
use xfw_konst_var
use xfw_netz_variablen
use xfw_functions

implicit none

! Variablendeklaration
character(len=*),                              intent(in)             :: input_file
logical,                                       intent(in),  optional  :: readmodel
logical,                                       intent(in),  optional  :: inviscid
double precision, dimension(:,:),              intent(in),  optional  :: displacements
double precision,                              intent(in),  optional  :: dispfac
double precision,                              intent(in),  optional  :: given_dalpha
logical,                                       intent(in),  optional  :: paneling
double precision, dimension(:,:), allocatable, intent(out), optional  :: pressures
double precision, dimension(:), allocatable,   intent(out), optional  :: F_total
double precision, dimension(:), allocatable,   intent(out), optional  :: al
double precision, dimension(:), allocatable,   intent(out), optional  :: ca
double precision, dimension(:), allocatable,   intent(out), optional  :: cd
double precision, dimension(:), allocatable,   intent(out), optional  :: cm
double precision, dimension(:), allocatable,   intent(out), optional  :: hkmax
double precision, dimension(:), allocatable,   intent(out), optional  :: Macrit
logical,                                       intent(in),  optional  :: output

double precision                                                      :: C_lift,C_drag,C_mom
double precision                                                      :: H_max,mach_crit
double precision                                                      :: helpfac
double precision                                                      :: xfoil_val
double precision                                                      :: start_alpha
double precision                                                      :: da_restart = 2.5d-1
double precision, dimension(:), allocatable                           :: cp,p_mano

logical                                                               :: xf_visc
logical                                                               :: success
logical                                                               :: xf_plot=.false.
logical                                                               :: xf_pane=.false.
logical                                                               :: force_inviscid=.false.
logical                                                               :: read_model=.true.
logical                                                               :: write_output=.true.

character(len=7)                                                      :: prefix_af = 'xf_arfl'
character(len=7)                                                      :: prefix_rs = 'xf_rslt'
character(len=7)                                                      :: prefix_hk = 'xf_rshk'
character(len=7)                                                      :: prefix_cp = 'xf_rscp'
character(len=7)                                                      :: prefix_xf = 'xf_inpt'

character(len=256)                                                    :: filename_af
character(len=256)                                                    :: filename_rs
character(len=256)                                                    :: filename_hk
character(len=256)                                                    :: filename_cp
character(len=256)                                                    :: filename_xf

integer                                                               :: ii
integer                                                               :: n_call=0
integer                                                               :: xf_iter=1000
integer                                                               :: restart_max=15

interface

    subroutine xfw_write_xfoil_input(plot,re,ma,n_crit,xtr_up,xtr_lo,iterXfoil,visc,pane,autotrim,value,fname_af,fname_rs,fname_hk,fname_cp,fname,start_alpha)
        implicit none
        logical                , intent(in)           :: plot
        double precision       , intent(in)           :: re
        double precision       , intent(in)           :: ma
        integer                , intent(in)           :: n_crit
        double precision       , intent(in)           :: xtr_up
        double precision       , intent(in)           :: xtr_lo
        integer                , intent(in)           :: iterXfoil
        logical                , intent(in)           :: visc
        logical                , intent(in)           :: pane
        logical                , intent(in)           :: autotrim
        double precision       , intent(in)           :: value
        character(len=256)     , intent(in)           :: fname_af
        character(len=256)     , intent(in)           :: fname_rs
        character(len=256)     , intent(in)           :: fname_hk
        character(len=256)     , intent(in)           :: fname_cp
        character(len=256)     , intent(in)           :: fname
        double precision       , intent(in), optional :: start_alpha
    end subroutine xfw_write_xfoil_input
    
    subroutine xfw_read_xfoil(alpha,CL,CD,CM,max_HK,cpx,mach_crit,success,pane,visc,fname_rs,fname_hk,fname_cp)
        use xfw_netz_variablen
        implicit none
        double precision       , intent(out)                        :: alpha
        double precision       , intent(out)                        :: CL
        double precision       , intent(out)                        :: CD
        double precision       , intent(out)                        :: CM
        double precision       , intent(out)                        :: max_HK
        double precision       , intent(out), dimension(n_nodes)    :: cpx
        double precision       , intent(out)                        :: mach_crit
        logical                , intent(out)                        :: success
        logical                , intent(in)                         :: pane
        logical                , intent(in)                         :: visc
        character(len=256)     , intent(in)                         :: fname_rs
        character(len=256)     , intent(in)                         :: fname_hk
        character(len=256)     , intent(in)                         :: fname_cp
    end subroutine

end interface

! Beginne Programmausfuehrung

if (present(readmodel)) then
    read_model = readmodel
end if

if (present(inviscid)) then
    force_inviscid = inviscid
end if

if (present(output)) then
    write_output = output
end if

if (present(given_dalpha)) then
    aoa = aoa + given_dalpha
    cl_given = .false.
end if

if (present(paneling)) then
    xf_pane = paneling
end if

! Lese Eingabedateien
if (read_model .EQV. .TRUE.) then
    if (allocated(nodes)) deallocate(nodes)

    call xfw_read_input(trim(input_file) // '.dat')
end if

! do not read model again on next call of subroutine
read_model = .false.

! add displacements to nodes
if (present(displacements)) then
    helpfac = 1.d0/chordscale
    if (present(dispfac)) then
        helpfac = helpfac*dispfac
    end if
    if (size(displacements,1) .ne. (n_nodes - 1)) then
        write(*,*) 'Wrong number of displacements given!'
        stop
    end if
    if (invert_points .eqv. .false.) then
        do ii = 1, size(displacements,1)
            nodes(ii,1) = nodes(ii,1) + displacements(ii,1)*helpfac
            nodes(ii,2) = nodes(ii,2) + displacements(ii,2)*helpfac
        end do
        nodes(n_nodes,1) = nodes(n_nodes,1) + displacements(1,1)*helpfac
        nodes(n_nodes,2) = nodes(n_nodes,2) + displacements(1,2)*helpfac
    else
        do ii = 2, size(displacements,1)+1
            nodes(ii,1) = nodes(ii,1) + displacements(size(displacements,1)-ii+2,1)*helpfac
            nodes(ii,2) = nodes(ii,2) + displacements(size(displacements,1)-ii+2,2)*helpfac
        end do
        nodes(1,1) = nodes(1,1) + displacements(1,1)*helpfac
        nodes(1,2) = nodes(1,2) + displacements(1,2)*helpfac
    end if
end if

! calculate panel lengths
allocate (lengths(n_panels))
do ii=1,n_panels
    lengths(ii) = abs_vector2(nodes(ii+1,:)-nodes(ii,:))
end do

if ((reynolds .gt. 0.d0) .and. (force_inviscid .eqv. .false.)) then
    xf_visc = .true.
else
    xf_visc = .false.
end if

if (cl_given .eqv. .true.) then
    xfoil_val = cl_target
else
    xfoil_val = aoa
end if

n_call = n_call + 1

! write airfoil coordinates
write(filename_af,'(A<len_trim(prefix_af)>,A,I0.5,A4)') TRIM(prefix_af),  '_', n_call, '.dat'
call xfw_write_af_coords(filename_af)

! allocate arrays
allocate (cp(n_nodes),p_mano(n_nodes))

! Wenn Xfoil nicht konvergierte, Neustart mit leicht variiertem Anstellwinkel als Startwert der Iteration
do ii = 1, restart_max

    write(*,*) '    XFOIL-START: ', ii, '/', restart_max
    
    ! write xfoil-input
    write(filename_rs,'(A<len_trim(prefix_rs)>,A,I0.5,A4)') TRIM(prefix_rs), '_', n_call, '.dat'
    write(filename_hk,'(A<len_trim(prefix_hk)>,A,I0.5,A4)') TRIM(prefix_hk), '_', n_call, '.dat'
    write(filename_cp,'(A<len_trim(prefix_cp)>,A,I0.5,A4)') TRIM(prefix_cp), '_', n_call, '.dat'
    write(filename_xf,'(A<len_trim(prefix_xf)>,A,I0.5,A4)') TRIM(prefix_xf), '_', n_call, '.dat'

    if (ii .eq. 1) then
        call xfw_write_xfoil_input(xf_plot,reynolds,Ma,n_crit,xtr_up,xtr_lo,xf_iter,xf_visc,xf_pane,cl_given,xfoil_val,filename_af,filename_rs,filename_hk,filename_cp,filename_xf)
    else
        if (cl_given .eqv. .true.) then
            ! Leichte Aenderung des Anstellwinkels als Startwert in Xfoil, wenn vorherige Rechnung nicht konvergiert ist
            start_alpha = dble(ii-1)*da_restart
        else
            ! Leichte Aenderung des Anstellwinkels als Startwert in Xfoil, wenn vorherige Rechnung nicht konvergiert ist
            start_alpha = xfoil_val - dble(ii-1)*da_restart
        end if

        call xfw_write_xfoil_input(xf_plot,reynolds,Ma,n_crit,xtr_up,xtr_lo,xf_iter,xf_visc,xf_pane,cl_given,xfoil_val,filename_af,filename_rs,filename_hk,filename_cp,filename_xf,start_alpha)
    end if

    call xfw_start_xfoil(filename_xf)

    call xfw_read_xfoil(aoa,C_lift,C_drag,C_mom,H_max,cp,mach_crit,success,xf_pane,xf_visc,filename_rs,filename_hk,filename_cp)

    ! Loesche Dateien
    call SYSTEM('rm ' // TRIM(filename_rs) // ' ' // TRIM(filename_cp) // ' ' // TRIM(filename_xf))
    if (xf_visc .eqv. .true.) then
      call SYSTEM('rm ' // TRIM(filename_hk))
    end if
    
    if (success .eqv. .true.) exit
    
!     if (cl_given .eqv. .true.) exit
    
    n_call = n_call + 1
end do
        
if (success .eqv. .false.) then
    write(*,*) 'Error in Xfoil-calculation!'
    STOP
end if

! Loesche Profildatei
call SYSTEM('rm ' // TRIM(filename_af))

! return manometer pressure values
if (present(pressures) .and. (xf_pane .eqv. .true.)) then
    write(*,*) 'Rueckgabe von Druckwerten bei Neuvernetzung in Xfoil nicht moeglich!'
    STOP
end if
p_mano(:) = cp(:)*q_infty
if (present(pressures)) then
    allocate(pressures(n_panels,1))
    if (invert_points .eqv. .false.) then
        do ii = 1, n_panels
             pressures(ii,1) = (p_mano(ii) + p_mano(ii+1))/2.d0
        end do
    else
        do ii = 1, n_panels
            pressures(ii,1) = (p_mano(1 + n_panels - ii) + p_mano(2 + n_panels - ii))/2.d0
        end do
    end if
end if

! return total aerodynamical force
if (present(F_total)) then
    allocate(F_total(1))
    F_total(1) = 0.d0
    do ii = 1, n_panels
        F_total(1) = F_total(1) + abs(pressures(ii,1))*lengths(ii)
    end do
    F_total(1) = F_total(1)*chordscale
end if

! return critical mach number
if (present(Macrit)) then
    allocate(Macrit(1))
    Macrit(1) = mach_crit
end if

! return angle of attack in degrees
if (present(al)) then
    allocate(al(1))
    al = aoa
end if

! return lift coefficient
if (present(ca)) then
    allocate(ca(1))
    ca = C_lift
end if

! return drag coefficient
if (present(cd)) then
    allocate(cd(1))
    cd = C_drag
end if

! return moment coefficient
if (present(cm)) then
    allocate(cm(1))
    cm = C_mom
end if

! return maximum kinematic shape factor
if (present(hkmax)) then
    allocate(hkmax(1))
    hkmax = H_max
end if

! remove displacements from nodes
if (present(displacements)) then
    helpfac = 1.d0/chordscale
    if (present(dispfac)) then
        helpfac = helpfac*dispfac
    end if
    if (invert_points .eqv. .false.) then
        do ii = 1, size(displacements,1)
            nodes(ii,1) = nodes(ii,1) - displacements(ii,1)*helpfac
            nodes(ii,2) = nodes(ii,2) - displacements(ii,2)*helpfac
        end do
        nodes(n_nodes,1) = nodes(n_nodes,1) - displacements(1,1)*helpfac
        nodes(n_nodes,2) = nodes(n_nodes,2) - displacements(1,2)*helpfac
    else
        do ii = 2, size(displacements,1)+1
            nodes(ii,1) = nodes(ii,1) - displacements(size(displacements,1)-ii+2,1)*helpfac
            nodes(ii,2) = nodes(ii,2) - displacements(size(displacements,1)-ii+2,2)*helpfac
        end do
        nodes(1,1) = nodes(1,1) - displacements(1,1)*helpfac
        nodes(1,2) = nodes(1,2) - displacements(1,2)*helpfac
    end if
end if

! Deallocate Arrays
if (allocated(lengths)) deallocate(lengths)
if (allocated(cp))      deallocate(cp)
if (allocated(p_mano))  deallocate(p_mano)

end subroutine xfoilwrapper_main
