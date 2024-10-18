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
subroutine panel2d_main(input_file, readmodel, inviscid, displacements, dispfac, given_dalpha, paneling, pressures, F_total, al, ca, cd, cm, hkmax, xsep, Mcrit, output)

! Load modules
use konst_var
use netz_variablen
use functions

implicit none

interface
    subroutine add_displacements(displacements, scle)
        double precision, dimension(:,:), intent(in)  :: displacements
        double precision,                 intent(in)  :: scle
    end subroutine add_displacements
end interface
 
interface
    subroutine pressure_calculation(ue, cp, p_mano, C_lift, C_mom, Ma_crit)
        use netz_variablen
        double precision, dimension(n_panels), intent(in)  :: ue
        double precision, dimension(n_panels), intent(out) :: cp, p_mano
        double precision, intent(out)                      :: C_lift, C_mom
        double precision, intent(out), optional            :: Ma_crit
    end subroutine pressure_calculation
end interface

! Declaration of variables
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
double precision, dimension(:), allocatable,   intent(out), optional  :: xsep
double precision, dimension(:), allocatable,   intent(out), optional  :: Mcrit
logical,                                       intent(in),  optional  :: output

double precision, dimension(2)                                        :: uw_infty,cl_temp
double precision, dimension(3)                                        :: aoa_try
double precision, dimension(:), allocatable                           :: bvec,sigma,sigma_invisc,cp,p_mano,ue,uen,phi,sigma_vii
type(bl_point), dimension(:), allocatable                             :: bl_values
double precision, dimension(:,:), allocatable                         :: Amat,Bmat,Amat_tmp
double precision                                                      :: C_lift,C_mom,C_drag,DeltaCL,FactorCM,Ma_crit
double precision                                                      :: dxsep,H_max,temp,transloc_up,transloc_lo
double precision, dimension(:,:), allocatable                         :: nodes_save
double precision                                                      :: cl_eps=1.d-4,cl_damp=1.d0
double precision                                                      :: helpfac
double precision                                                      :: af_scale,af_rotation

logical                                                               :: verbose        = .TRUE.
logical                                                               :: remesh         = .FALSE.
logical                                                               :: force_inviscid = .FALSE.
logical                                                               :: read_model     = .TRUE.
logical                                                               :: write_output   = .FALSE.
logical                                                               :: trim_cl
logical                                                               :: bl_error

integer                                                               :: ii
integer                                                               :: iter_cl
integer                                                               :: max_iter_cl=1000
integer                                                               :: info
integer, dimension(:), allocatable                                    :: ipiv

!*****************************************************************************!
! PANEL2D - 2D Panel program for the aerodynamic analysis of airfoils         !
!                                                                             !
! The theory of the inviscid panel program is based on:                       !
! [1]: p. 290ff (11.3.1 Combined Source and Doublet Method) and the           !
! corresponding example program D9 [1, p. 567ff].                             !
!                                                                             !
! The boundary layer calculation is based on:                                 !
! [2]: p. 191ff (7 The Boundary Layer) and the corresponding example program  !
! INTGRL [2, p. 255ff].                                                       !
!                                                                             !
! [1] Katz, J.; Plotkin, A.: Low-Speed Aerodynamics, 2nd ed., Cambridge:      !
!     Cambridge University Press, 2001, ISBN: 0-521-66552-3                   !
! [2] Moran, J.: An Introduction to Theoretical and Computational             !
!     Aerodynamics, New York: John Wiley & Sons, 1984, ISBN: 0-471-87491-4    !
!                                                                             !
! (c) 2018-2020, Florian Dexl                                                 !
!*****************************************************************************!


! Begin program execution

if (present(readmodel))   read_model     = readmodel
if (present(inviscid))    force_inviscid = inviscid
if (present(paneling))    remesh         = paneling
if (present(output))      write_output   = output

if (present(given_dalpha)) then
    aoa = aoa + deg_to_rad(given_dalpha)
    cl_given = .false.
end if

! Read input file
if (read_model .EQV. .TRUE.) then
    if (allocated(nodes)) deallocate(nodes)
    call read_input(TRIM(input_file) // '.dat', verbose)

    ! Do not read model again on next call of subroutine
    read_model = .FALSE.
end if

! Name of output files as input file
res_prefix = TRIM(input_file) // '_'

! Save original node coordinates
allocate(nodes_save(n_nodes,2))
nodes_save(:,:) = nodes(:,:)

! Add displacements to nodes
if (present(displacements)) then
    helpfac = 1.d0/chordscale
    if (present(dispfac)) then
        helpfac = helpfac*dispfac
    end if
    call add_displacements(displacements, helpfac)
end if

! Check input data
call check_airfoil()

! Unification of airfoil
call unify_airfoil(af_scale, af_rotation)
! Adapt given angle of attack to rotation
aoa = aoa - af_rotation
! Adapt Reynolds number
reynolds = reynolds/af_scale

! Remesh using Xfoil's meshing routine
if (remesh .EQV. .TRUE.) then
    call refine_mesh(verbose)
end if

! Calculate geometric values
call calculate_geometry()

! Calculate panel geometry
call panel_geometry()

! Allocate arrays
allocate(Amat(n_nodes,n_nodes),  &
       & Bmat(n_panels,n_panels) )
allocate(sigma(n_panels), sigma_invisc(n_panels), phi(n_panels), sigma_vii(n_panels))
allocate(cp(n_panels), ue(n_panels), p_mano(n_panels))
allocate(Amat_tmp(n_nodes,n_nodes))
allocate(ipiv(n_nodes))
allocate(bvec(n_nodes))
allocate(uen(n_nodes))
allocate(bl_values(n_nodes))

bl_values(:)%H      = 0.d0
bl_values(:)%theta  = 0.d0
bl_values(:)%deltas = 0.d0
bl_values(:)%cf     = 0.d0
bl_values(:)%nfac   = 0.d0

! Initial angles of attack for trimming to target-C_L value
aoa_try(1) = deg_to_rad(0.d0)  - af_rotation
aoa_try(2) = deg_to_rad(0.5d0) - af_rotation

! Calculate A and B matrices
call calculate_ABmat(Amat, Bmat)

! Perform LU factorization of A matrix
call dgetrf(n_nodes, n_nodes, Amat, n_nodes, ipiv, info)
if (info .ne. 0) then
    write(*,*) 'ERROR in DGETRF. Returned info: ', info
    stop
end if
! In Amat is now the LU factorization of the A matrix

! Loop for trimming to target-C_L value
iter_cl = 0
trim_cl = cl_given
do while (.true.)
    iter_cl = iter_cl + 1

    if (trim_cl .eqv. .true.) then
        if (iter_cl .le. 2) then
            aoa = aoa_try(iter_cl)
        else
            aoa = aoa_try(3)
        end if
    end if

    ! Calculate velocity vector of outer flow
    uw_infty(1) = cos(aoa)
    uw_infty(2) = sin(aoa)

    ! Calculate source strengths
    do ii = 1, n_panels
        sigma_invisc(ii) = DOT_PRODUCT(uw_infty(:),-1.d0*n_vec(ii,:))
    end do

    sigma_vii(:)     = 0.d0

    ! Add boundary condition of viscous solution to source strengths
    sigma(:) = sigma_invisc(:) - sigma_vii(:)

    ! Calculate RHS vector
    bvec(1:n_panels) = matmul(Bmat(:,:),sigma(:))
    bvec(n_panels+1) = 0.d0

    ! Solve for doublet strengths
    call dgetrs('N', n_nodes, 1, Amat, n_nodes, ipiv, bvec, n_nodes, info)
    if (info .ne. 0) then
        write(*,*) 'ERROR in DGETRS. Returned info: ', info
        stop
    end if
    ! In bvec are now doublet strengths

    ! Calculcate aerodynamic results
    do ii = 1, n_panels
        phi(ii) = DOT_PRODUCT(coll(ii,:),uw_infty(:)) + bvec(ii)
    end do

    ! Calculate tangential velocity at nodes
    do ii = 1, n_panels-1
        uen(ii+1) = -1.d0*(phi(ii) - phi(ii+1))/((lengths(ii+1) + lengths(ii))/2.d0)
    end do

    ! Still the velocities at the trailing edge are missing.
    ! According to [2, S. 94], for fulfilling the Kutta condition,
    ! the velocities at the upper and lower side of the airfoil at
    ! the trailing edge must have an identical absolute value.
    ! Here, the velocities at the trailing edge are extrapolated first
    ! and the absolute value is averaged from the upper and lower side.
    ! The sign of the upper and lower side is kept, thereby.
    uen(1) = (uen(2) - uen(3))/lengths(2)*(lengths(2) + lengths(1)) + uen(3)
    uen(n_nodes) = (uen(n_nodes-1) - uen(n_nodes-2))/lengths(n_panels-1)*(lengths(n_panels-1) + lengths(n_panels)) + uen(n_nodes-2)

    ! Interpolate velocity at collocation points
    do ii = 1, n_panels
        ue(ii) = (uen(ii) + uen(ii+1))/2.d0
    end do

    ! Pressure calculation
    if (present(Mcrit)) then
      call pressure_calculation(ue,cp,p_mano,C_lift,C_mom,Ma_crit)
    else
      call pressure_calculation(ue,cp,p_mano,C_lift,C_mom)
      Ma_crit = 0.d0
    end if
    
    ! Boundary layer calculation
    if ((reynolds .gt. 0.d0) .and. (force_inviscid .eqv. .false.)) then
        call boundary_layer(verbose,n_nodes,nodes,uen,reynolds,C_drag,DeltaCL,FactorCM,dxsep,transloc_up,transloc_lo,bl_values,bl_error,C_lift)
    else
        C_drag   = 0.d0
        DeltaCL  = 0.d0
        FactorCM = 1.d0
        dxsep    = 0.d0
        transloc_up = 0.d0
        transloc_lo = 0.d0
        bl_error = .false.
    end if

    if (bl_error .eqv. .true.) then
        write(*,*) 'Error in boundary layer calculation.'
        stop
    end if

    C_lift = C_lift*(1.d0 + DeltaCL)

    C_mom  = C_mom*FactorCM

    if (trim_cl .eqv. .false.) then
        exit
    else
        if (iter_cl .le. 2) then
            cl_temp(iter_cl) = C_lift
        else if (iter_cl .gt. max_iter_cl) then
            write(*,*)
            write(*,*) 'Iteration alpha(cl) was not successful.'
            STOP
        else
            if (abs(C_lift - cl_target) .lt. cl_eps) then
                aoa_try(3) = aoa
                trim_cl    = .false.
                write(*,*)
                write(*,*) 'Alpha found for target lift coefficient: ', (aoa*180.d0/pi), ' degrees.'
                write(*,'(A12,F6.2,A10)') 'CL_EPS ONLY ', (cl_eps*1.d2), ' PERCENT !'
                exit
            end if
            cl_temp(1) = cl_temp(2)
            cl_temp(2) = C_lift
            aoa_try(1) = aoa_try(2)
            aoa_try(2) = aoa_try(3)
        end if
        if (iter_cl .ge. 2) then
            aoa_try(3) = aoa_try(2) - (cl_temp(2) - cl_target)*(aoa_try(2) - aoa_try(1))/(cl_temp(2) - cl_temp(1))*cl_damp
        end if
    end if
end do

if (dxsep .gt. 0.d0) then
   write(*,*) 'Separation occured, dsxep = ', dxsep
end if

!if (dxsep .gt. 0.01d0) then
!    write(*,*) "Separation over more than 1% of the airfoil's surface"
!    write(*,*) 'Non trustable values at final calculation!'
!    write(*,*) 'Penalty drag coefficient!'
!    C_drag = 1.d5
!end if

! Return manometer pressure values
if (present(pressures)) then
    if (remesh .EQV. .TRUE.) then
        write(*,*) 'ERROR: Pressures cannot be returned with remeshing enabled!'
        STOP
    else
        allocate(pressures(n_panels,1))
        if (invert_points .eqv. .false.) then
            pressures(:,1) = p_mano(:)
        else
            pressures(:,1) = p_mano(n_panels:1:-1)
        end if
    end if
end if

if (present(F_total)) then
    allocate(F_total(1))
    F_total(1) = sum(abs(p_mano(:))*lengths(:))
    F_total(1) = F_total(1)*chordscale
end if

! Reset angle of attack by rotation angle
aoa = aoa + af_rotation

! Return angle of attack in degrees
if (present(al)) then
    allocate(al(1))
    al = rad_to_deg(aoa)
end if

! Return lift coefficient
if (present(ca)) then
    allocate(ca(1))
    if (C_lift /= C_lift) then
        ca = 1.d7
    else
        ca = C_lift
    end if
end if

! Return drag coefficient
if (present(cd)) then
    allocate(cd(1))
    if (C_drag /= C_drag) then
        cd = 1.d7
    else
        cd = C_drag
    end if
end if

! Return moment coefficient
if (present(cm)) then
    allocate(cm(1))
    if (C_mom /= C_mom) then
        cm = 1.d7
    else
        cm = C_mom
    end if
end if

! Return maximum kinematic shape factor
H_max = maxval(bl_values(:)%H)
if (present(hkmax)) then
    allocate(hkmax(1))
    if (H_max /= H_max) then
        hkmax = 1.d7
    else
        hkmax = H_max
    end if
end if

! Return separation position
if (present(xsep)) then
    allocate(xsep(1))
    xsep = 1.d0 - dxsep
end if

! Return critical Mach-Number
if (present(Mcrit)) then
    allocate(Mcrit(1))
    Mcrit = Ma_crit
end if


if (write_output .eqv. .true.) then
    call write_vtk(TRIM(res_prefix)//'out.vtk',chordscale,cp,ue,uen,p_mano,bl_values)
end if

! Restore original node coordinates
deallocate(nodes)
n_nodes = size(nodes_save,1)
allocate(nodes(n_nodes,2))
nodes(:,:) = nodes_save(:,:)
deallocate(nodes_save)
! Restore original Reynolds number
reynolds = reynolds*af_scale

! Wenn C_drag ungleich NaN ist, Ergebnisdatei schreiben
if ((.not. (C_drag /= C_drag)) .and. (bl_error .eqv. .false.) .and. (write_output .eqv. .true.)) then
    open(unit=21,file=TRIM(res_prefix)//TRIM(res_file),status='replace')
    write(21,'(A17)') '# Results Panel2D'
    write(21,'(E24.16,X,A13)') aoa*180.d0/pi, 'aoa [degrees]'
    write(21,'(E24.16,X,A7)')  C_lift,        'C_L [-]'
    write(21,'(E24.16,X,A7)')  C_drag,        'C_D [-]'
    write(21,'(E24.16,X,A7)')  C_mom,         'C_M [-]'
    write(21,'(E24.16,X,A9)')  H_max,         'H_max [-]'
    write(21,'(E24.16,X,A11)') Ma_crit,       'Ma_crit [-]'
    write(21,'(E24.16,X,A9)')  dxsep,         'dxsep [-]'
    close(21)
end if

if (verbose .eqv. .true.) then
    write(*,*)
    write(*,'(A7,8X,E24.16)')  'C_lift:',   C_lift
    write(*,'(A7,8X,E24.16)')  'C_drag:',   C_drag
    write(*,'(A6,9X,E24.16)')  'C_mom:',    C_mom
    write(*,'(A6,9X,E24.16)')  'H_max:',    H_max
    write(*,'(A4,11X,E24.16)') 'L/D:',      C_lift/C_drag
    write(*,'(A8,7X,E24.16)')  'Ma_crit:',  Ma_crit
    write(*,'(A6,9X,E24.16)')  'dxsep:',    dxsep
    write(*,'(A11,4X,E24.16)')  'transloc_up:', transloc_up/100.d0
    write(*,'(A11,4X,E24.16)')  'transloc_lo:', transloc_lo/100.d0
end if

! Deallokiere Arrays
if (allocated(Amat))          deallocate(Amat)
if (allocated(ipiv))          deallocate(ipiv)
if (allocated(Bmat))          deallocate(Bmat)
if (allocated(bvec))          deallocate(bvec)
if (allocated(sigma))         deallocate(sigma)
if (allocated(sigma_invisc))  deallocate(sigma_invisc)
if (allocated(sigma_vii))     deallocate(sigma_vii)
if (allocated(phi))           deallocate(phi)
if (allocated(cp))            deallocate(cp)
if (allocated(ue))            deallocate(ue)
if (allocated(n_vec))         deallocate(n_vec)
if (allocated(t_vec))         deallocate(t_vec)
if (allocated(lengths))       deallocate(lengths)
if (allocated(alpha_i))       deallocate(alpha_i)
if (allocated(coll))          deallocate(coll)
if (allocated(bl_values))     deallocate(bl_values)
if (allocated(uen))           deallocate(uen)

end subroutine panel2d_main
