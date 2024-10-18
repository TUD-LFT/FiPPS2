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
! Copyright (C) 2010-2012  Daniel Filkovic

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %                                                      %
! %            APAME - Aircraft Panel Method             %
! %______________________________________________________%
! %                                                      %
! %              3D potential flow solver                %
! %                                                      %
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! This file is part of APAME.

! APAME is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! APAME is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with APAME.  If not, see <http://www.gnu.org/licenses/>.

! file subr_apame_main.f90

subroutine apame_main(fname, readmodel, displacements, dispfac, pressures, F_total, cdi, output)

use module_kind_and_konst
use module_input
use module_grid
use module_modify_wake
use module_influence
use module_rhs
use module_velocity
use module_pressure
use module_induced_drag

implicit none

  ! INTENTS IN ===================================================================
  character(len=*), optional                                                 :: fname
  logical,                                            intent(in),  optional  :: readmodel
  real(kind=kind_float), dimension(:,:),              intent(in),  optional  :: displacements
  real(kind=kind_float),                              intent(in),  optional  :: dispfac
  real(kind=kind_float), dimension(:,:), allocatable, intent(out), optional  :: pressures
  real(kind=kind_float), dimension(:), allocatable,   intent(out), optional  :: F_total
  real(kind=kind_float), dimension(:), allocatable,   intent(out), optional  :: cdi
  logical,                                            intent(in),  optional  :: output
  ! PRIVATE ======================================================================
  real(kind=kind_float)                                                      :: helpfac
  integer                                                                    :: err
  integer                                                                    :: ii
  integer                                            :: interactive                    ,&
                                                      & info_solver                    ,&
                                                      & slash_position                 ,&
                                                      & alloc_stat
  real(kind=kind_float)                              :: start_time
  integer,      allocatable, dimension(:,:)          :: ipiv
  real(kind=kind_float), allocatable, dimension(:)   :: speed_x,speed_y,speed_z
  character(len=100)                                 :: file_name, job_name
  character(len=3)                                   :: version, subversion
  character(len=6)                                   :: build_date
  logical                                            :: read_model = .true.
  logical                                            :: file_exists
  logical                                            :: apame_error
  
  real(kind=kind_float), dimension(3)                :: aoa_try
  
  character(256)                                     :: message
  
  character(100)                                     :: lastfname
  
  parameter (version="3")             ! program version
  parameter (subversion="1")          ! program subversion
  parameter (build_date="180404_lft") ! build date
  ! ==============================================================================

  apame_error = .TRUE.

  if (present(readmodel)) then
    read_model = readmodel
  end if
  
  ! READING APAME INPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  if (read_model .eqv. .true.) then
        if (allocated(panel_type) .eq. .true.) deallocate(panel_type)
        if (allocated(node1) .eq. .true.) deallocate(node1)
        if (allocated(node2) .eq. .true.) deallocate(node2)
        if (allocated(node3) .eq. .true.) deallocate(node3)
        if (allocated(node4) .eq. .true.) deallocate(node4)
        if (allocated(elem1) .eq. .true.) deallocate(elem1)
        if (allocated(elem2) .eq. .true.) deallocate(elem2)
        if (allocated(elem3) .eq. .true.) deallocate(elem3)
        if (allocated(elem4) .eq. .true.) deallocate(elem4)
        if (allocated(alfa) .eq. .true.) deallocate(alfa)
        if (allocated(beta) .eq. .true.) deallocate(beta)
        if (allocated(x_node) .eq. .true.) deallocate(x_node)
        if (allocated(y_node) .eq. .true.) deallocate(y_node)
        if (allocated(z_node) .eq. .true.) deallocate(z_node)
  
      if (present(fname)) then
        file_name = fname
        interactive=0
      else
        ! getting command line argument
        call get_command_argument(1,file_name)
        ! in case no input file is given as argument, APAME will ask for one and enter interactive mode
        if (file_name .eq. "") then
            ! entering interactive mode
            interactive=1

            ! printing title
            call func_title( version                                                ,&
                            & subversion                                            ,&
                            & build_date                                             )
            
            ! requesting user for file name
            read *, file_name
        else
            ! entering non-interactive mode CHANGE THIS VALUE TO "0" WHEN RELEASING
            interactive=0
        endif
      endif
  else
      if (present(fname)) then
        file_name = fname
        interactive=0
      endif  
  endif

  ! create (replace) log file
  
      
  if (present(output)) then
    if (output .eq. .false.) then
      res_req = 0
      vtk_req = 0
      open(unit=2,file="/dev/null")
    else
      open(unit=2,file=trim(file_name)//".log",status="replace")
    end if
  else
      open(unit=2,file=trim(file_name)//".log",status="replace")
  endif

  if (read_model .eqv. .true.) then

      call func_tic(start_time)
      call func_new_line(interactive)

      ! check if input file exists
      inquire(file=trim(file_name)//".inp",exist=file_exists)
  
      ! open input file for usage
      if (file_exists .neqv. .true.) then
          call func_message( interactive                                          ,&
                           & "    ERROR: Input file "//trim(file_name)//".inp not found")
          goto 999
      endif
  
      ! print date and time
      call subr_write_time(interactive)
      call func_new_line(interactive)
  
      ! reading input file
      call input(interactive, file_name, version, subversion)
      if (input_err .ne. 0) then
          goto 999
      endif
      
      if (present(output)) then
        if (output .eq. .false.) then
          res_req = 0
          vtk_req = 0
        end if
      endif 
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      
      ! print jobname
      job_name=trim(file_name)
      ! find slash or backslash position if file was given with path
      slash_position=scan(job_name, '\/', back=.true.)
      if (slash_position .eq. 0) then
          call func_message( interactive                                          ,&
                           & "    Job name: "//trim(job_name))
      else
          call func_message( interactive                                          ,&
                           & "    Job name: "//job_name(slash_position+1:len_trim(job_name)))
      endif
      
      call func_new_line(interactive)
      call func_new_line(interactive)
        
      call func_toc( interactive, start_time                                      ,&
                   & "    Time for reading input.............................:")
      call func_new_line(interactive)
  end if

  ! do not read model again on next call of subroutine
  read_model = .false.
  
  ! add displacements to nodes
  if (present(displacements)) then
      if (present(dispfac)) then
          helpfac = dispfac
      else
          helpfac = 1._kind_float
      end if
      do ii = 1, size(displacements,1)
          x_node(ii) = x_node(ii) + displacements(ii,1)*helpfac
          y_node(ii) = y_node(ii) + displacements(ii,2)*helpfac
          z_node(ii) = z_node(ii) + displacements(ii,3)*helpfac
      end do
  end if
  
  ! print number of panels
  call subr_write_ptot(interactive, panel_num, panel_num_no_wake)
  
  ! print keywords
  call subr_write_key( interactive                                            ,&
                     & speed                                                  ,&
                     & ro                                                     ,&
                     & p_ref                                                  ,&
                     & mach                                                   ,&
                     & wing_span                                              ,&
                     & mac                                                    ,&
                     & wing_surf                                              ,&
                     & origin                                                 ,&
                     & method                                                 ,&
                     & symmetry                                               ,&
                     & wake_mod                                               ,&
                     & trefftz                                                ,&
                     & error                                                  ,&
                     & coplanarity_angle                                      ,&
                     & farfield                                               ,&
                     & collcalc                                               ,&
                     & velorder                                               )
  
  ! print cases (angles of attack and sideslip angles)
  call subr_write_case(interactive, case_num, alfa, beta, trim_cl, target_cl)
  call func_new_line(interactive)
  
  if (trim_cl .eq. 0) then
      ! Standard calculation
      call given_aoa(.false.,.true.,err)
      if (err .ne. 0) then
          goto 999
      end if
  else
      if (case_num .ne. 1) then
          call func_message( interactive                              ,&
                  & "ERROR: TRIM_CL only compatible with CASE_NUM = 1")
          goto 999
      end if
  
      ! Starting values for angle of attack 
      aoa_try(1) = -2._kind_float/180._kind_float*pi
      aoa_try(2) =  4._kind_float/180._kind_float*pi
      aoa_try(3) =  0._kind_float                         ! will be overwritten with final aoa
  
      ! Trim to target lift coefficient with given precision
      call given_cl(aoa_try, target_cl, cl_eps, .true., err)
      if (err .ne. 0) then
          goto 999
      end if
  end if

  ! return manometer pressure values
  if (present(pressures)) then
      allocate(pressures(panel_num_no_wake,case_num))
      pressures(:,:) = p_mano(:,:)
      write(*,*) 'target_cl        :', target_cl
      write(*,*) 'cl               :', coef_lift(1)
      write(*,*) 'cl_mom_corrected :', cl_mom_corrected(1)
  end if
  
  ! return summed force
  if (present(F_total)) then
      allocate(F_total(case_num))
      F_total(:) = Ftot(:)
  end if
  
  ! remove displacements from nodes
  if (present(displacements)) then
      if (present(dispfac)) then
          helpfac = dispfac
      else
          helpfac = 1._kind_float
      end if
      do ii = 1, size(displacements,1)
          x_node(ii) = x_node(ii) - displacements(ii,1)*helpfac
          y_node(ii) = y_node(ii) - displacements(ii,2)*helpfac
          z_node(ii) = z_node(ii) - displacements(ii,3)*helpfac
      end do
  end if
  
  if (present(cdi)) then
      allocate(cdi(case_num))
      cdi(:) = coef_idrag(:)
  end if

  apame_error = .FALSE.

999 continue

  if (apame_error .eqv. .TRUE.) then
       write(*,*) 'ERROR: APAME CALCULATION FAILED.'
       write(*,*) 'PLEASE READ APAME LOG-FILE ', trim(file_name), '.log .'
       STOP
  end if
      
  call func_new_line(interactive)
  
  close (2)
  
  ! deallocate arrays
  call deallocate_arrays()
  call deallocate_grid()

contains
  
  subroutine given_aoa(wake_only,write_res,err)
      !
      implicit none
      !
      ! INTENTS OUT===================================================================
      integer, intent(out)                                                       :: err
      logical, intent(in)                                                        :: wake_only,write_res
      ! ==============================================================================
      err = 0
      
      ! calculating airspeed components
      if (.not. allocated(speed_x)) allocate (speed_x(case_num))
      if (.not. allocated(speed_y)) allocate (speed_y(case_num))
      if (.not. allocated(speed_z)) allocate (speed_z(case_num))
      call subr_speeds( case_num                                                  ,&
                      & speed                                                     ,&
                      & alfa,beta                                                 ,&
                      & speed_x,speed_y,speed_z                                   )
      
      call calculate_mod_wake(alfa(1),beta(1),err)
      if (err .ne. 0) then
          goto 999
      end if
      
      call calculate_grid(wake_only,err)
      if (err .ne. 0) then
          goto 999
      end if
      
      call calculate_influence_coefficients(wake_only,err)
      if (err .ne. 0) then
          goto 999
      end if
      
      ! CALCULATING RIGHT-HAND SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      call func_tic(start_time)
      call func_message( interactive                                              ,&
                       & "Calculating right-hand side...")
      
      ! allocate b to minimum size for constant doublet method
      if ((method .eq. 1) .and. (.not. allocated(b))) allocate (b(1,1))
      
      call rhs_calc( b                                                            ,&
                   & method                                                       ,&
                   & panel_num                                                    ,&
                   & panel_num_no_wake                                            ,&
                   & interactive                                                  ,&
                   & case_num                                                     ,&
                   & speed_x,speed_y,speed_z                                      ,&
                   & n1,n2,n3                                                     ,&
                   & cx,cy,cz                                                     )
      if (rhs_err .ne. 0) then
          err = 1
          goto 999
      endif
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      
      call func_toc( interactive ,start_time                                      ,&
                   & "    Time for calculating right-hand side...............:")
      call func_new_line(interactive)
      
      ! SOLVING SYSTEM OF EQUATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      call func_tic(start_time)
      call func_message( interactive                                              ,&
                       & "Solving system of equations...")
      
      allocate ( ipiv(panel_num_no_wake, panel_num_no_wake)                       ,&
               & stat=alloc_stat                                                  )
      if (alloc_stat .ne. 0) then
          call func_message( interactive                                          ,&
                           & "    ERROR: Not enough memory for allocating permutation matrix &
                           & (main system of equations)")
          err = 1
          goto 999
      endif
      
      ! solving system of equations
      if (kind_float .eq. kind_float_l) then
          ! single precision
          call sgesv( panel_num_no_wake                                               ,&
                    & case_num                                                        ,&
                    & a                                                               ,&
                    & panel_num_no_wake                                               ,&
                    & ipiv                                                            ,&
                    & rhs                                                             ,&
                    & panel_num_no_wake                                               ,&
                    & info_solver                                                     )
      else
          ! double precision
          call dgesv( panel_num_no_wake                                               ,&
                    & case_num                                                        ,&
                    & a                                                               ,&
                    & panel_num_no_wake                                               ,&
                    & ipiv                                                            ,&
                    & rhs                                                             ,&
                    & panel_num_no_wake                                               ,&
                    & info_solver                                                     )
      end if
      
      deallocate(ipiv)
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      
      if (info_solver .lt. 0) then
          call func_message( interactive                                          ,&
                           & "    ERROR: Illegal value detected in system solver")
          err = 1
          goto 999
      elseif (info_solver .gt. 0) then
          call func_message( interactive                                          ,&
                           & "    ERROR: Factorization complete with singular solution")
          err = 1
          goto 999
      endif
      
      call func_new_line(interactive)
      call func_toc( interactive, start_time                                      ,&
                   & "    Time for solving system of equations...............:")
      call func_new_line(interactive)
      
      call calculate_postprocessing(case_num, alfa, beta, err)
      if (err .ne. 0) then
          goto 999
      end if

      if (write_res .eqv. .true.) then
          call calculate_results(err)
          if (err .ne. 0) then
              goto 999
          end if
      end if
      
999 continue
  
  end subroutine given_aoa


  subroutine given_cl(aoa, cl, eps, write_res, err)
      !
      implicit none
      !
      ! INTENTS OUT===================================================================
      integer, intent(out)                                                      :: err
      ! PRIVATE ======================================================================
      logical, intent(in)                                                       :: write_res
      real(kind=kind_float), dimension(2)                                       :: cl_tmp
      real(kind=kind_float), dimension(3), intent(inout)                        :: aoa
      real(kind=kind_float), intent(in)                                         :: cl, eps
      integer                                                                   :: alloc_stat,aoa_err
      integer                                                                   :: ii, iter, max_iter = 20
      logical                                                                   :: wake_only
      ! ==============================================================================
      err = 0
      !
      iter = 0

      do while (.true.)
          iter = iter + 1
      
          call func_new_line(interactive)
          write(message,'(A31,I0)') 'Iteration for finding aoa(cl): ', iter
          call func_message( interactive, trim(message))
          call func_new_line(interactive)
          call func_new_line(interactive)
      
          if (iter .le. 2) then
              alfa(:) = aoa(iter)    
          else
              alfa(:) = aoa(3)
          end if

          ! calculating airspeed components
          allocate ( speed_x(case_num)                                                ,&
                   & speed_y(case_num)                                                ,&
                   & speed_z(case_num)                                                )
          call subr_speeds( case_num                                                  ,&
                          & speed                                                     ,&
                          & alfa,beta                                                 ,&
                          & speed_x,speed_y,speed_z                                   )

          if (iter .eq. 1) then
              wake_only = .false.
          else
              wake_only = .true.
          end if
      
          call given_aoa(wake_only,.false.,aoa_err)
          if (aoa_err .ne. 0) then
              err = 1
              goto 999
          endif

          if (iter .le. 2) then
              cl_tmp(iter)  = cl_mom_corrected(1)
          else if (iter .gt. max_iter) then
              call func_message( interactive                                              ,&
                              & "Maximum number of iterations for finding aoa(cl) reached!")
              err = 1
              goto 999
          else
              if (abs(cl_mom_corrected(1) - cl) .lt. eps) then
                  exit
              end if
              cl_tmp(1)  = cl_tmp(2)
              cl_tmp(2)  = cl_mom_corrected(1)
              aoa(1) = aoa(2)
              aoa(2) = aoa(3)
          end if  
          aoa(3) = aoa(2) - (cl_tmp(2) - cl)*(aoa(2) - aoa(1))/(cl_tmp(2) - cl_tmp(1))

          ! deallocate arrays
          call deallocate_arrays()
      end do
      
      call func_new_line(interactive)
      write(message,'(A30,F7.4)') 'Target angle of attack [deg]: ', aoa_try(3)*180._kind_float/pi
      call func_message( interactive, trim(message))
      call func_new_line(interactive)
      
      if (write_res .eq. .true.) then
          call calculate_results(err)
          if (err .ne. 0) then
              goto 999
          end if
      end if

999 continue
  
  end subroutine given_cl

  
  subroutine calculate_grid(wake_only, err)
      !
      implicit none
      !
      ! INTENTS OUT ==================================================================
      logical, intent(in)                                                       :: wake_only
      integer, intent(out)                                                      :: err
      ! ==============================================================================
      err = 0
      !
      ! CALCULATING GRID INFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call func_tic(start_time)
      call func_message( interactive                                              ,&
                       & "Calculating grid information...")
      
      call grid( interactive                                                      ,&
               & collcalc                                                         ,&
               & panel_num                                                        ,&
               & panel_num_no_wake                                                ,&
               & node_num                                                         ,&
               & panel_type                                                       ,&
               & node1,node2,node3,node4                                          ,&
               & farfield                                                         ,&
               & error                                                            ,&
               & colldist                                                         ,&
               & x_node,y_node,z_node                                             ,&
               & wake_only                                                          )
      if (grid_err .ne. 0) then
          err = 1
          goto 999
      endif
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      
      call func_toc( interactive, start_time                                      ,&
                   & "    Time for preprocessing.............................:")
      call func_new_line(interactive)
      
999 continue

  end subroutine calculate_grid
  
  
  subroutine calculate_mod_wake(alfa_wake,beta_wake,err)
      !
      implicit none
      !
      ! INTENTS OUT ==================================================================
      real(kind=kind_float), intent(in)                                         :: alfa_wake, beta_wake
      integer, intent(out)                                                      :: err
      ! ==============================================================================
      err = 0
      !
      ! MODIFY WAKE PANELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call func_tic(start_time)
      call func_message( interactive                                              ,&
                       & "Modify wake panels.............")
      
      if (wake_mod .eq. 1) then
          if (case_num .ne. 1) then
              call func_message( interactive                              ,&
                               & "ERROR: WAKE_MOD only compatible with CASE_NUM = 1")
              err = 1
              goto 999
          else
              call modify_wake( interactive                                       ,&
                              & panel_num                                         ,&
                              & panel_num_no_wake                                 ,&
                              & node_num                                          ,&
                              & node1,node2,node3,node4                           ,&
                              & alfa_wake,beta_wake                               ,&
                              & trefftz(1), trefftz(2), trefftz(3)                ,&
                              & x_node,y_node,z_node                              )
          end if
      end if
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      
      call func_toc( interactive, start_time                                      ,&
                   & "    Time for modifying wake panels.....................:")
      call func_new_line(interactive)
  
999 continue
  
  end subroutine calculate_mod_wake


  subroutine calculate_influence_coefficients(wake_only, err)
      !
      implicit none
      !
      ! INTENTS OUT ==================================================================
      logical, intent(in)                                                       :: wake_only
      integer, intent(out)                                                      :: err
      ! ==============================================================================
      err = 0
      !
      ! CALCULATING INFLUENCE COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      call func_tic(start_time)
      call func_message( interactive                                              ,&
                       & "Calculating influence coefficients...")
      
      call influence( method                                                      ,&
                    & symmetry                                                    ,&
                    & panel_num                                                   ,&
                    & panel_num_no_wake                                           ,&
                    & interactive                                                 ,&
                    & panel_type                                                  ,&
                    & elem1,elem2                                                 ,&
                    & FF                                                          ,&
                    & S                                                           ,&
                    & error                                                       ,&
                    & colx,coly,colz                                              ,&
                    & cx,cy,cz                                                    ,&
                    & l1,l2,l3                                                    ,&
                    & p1,p2,p3                                                    ,&
                    & n1,n2,n3                                                    ,&
                    & x1,x2,x3,x4                                                 ,&
                    & y1,y2,y3,y4                                                 ,&
                    & d1,d2,d3,d4                                                 ,&
                    & wake_only                                                    )
      if (influence_err .ne. 0) then
          err = 1
          goto 999
      endif
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      
      call func_toc( interactive, start_time                                      ,&
                   & "    Time for calculating influence coefficients........:")
      call func_new_line(interactive)
  
999 continue
  
  end subroutine calculate_influence_coefficients

  
  subroutine calculate_postprocessing(case_num_tmp, alfa_tmp, beta_tmp, err)
      !
      implicit none
      !
      ! INTENTS OUT ==================================================================
      integer, intent(out)                                                      :: err
      ! PRIVATE ======================================================================
      integer, intent(in)                                                       :: case_num_tmp
      real(kind=kind_float), dimension(case_num_tmp), intent(in)                :: alfa_tmp, beta_tmp
      ! ==============================================================================
      err = 0
      !
      ! POSTPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
      call func_tic(start_time)
      call func_message( interactive                                              ,&
                      & "Postprocessing...")
      
      ! calculating induced velocities
      call velocities( rhs                                                        ,&
                      & node_num                                                   ,&
                      & panel_num                                                  ,&
                      & panel_num_no_wake                                          ,&
                      & case_num_tmp                                               ,&
                      & velorder                                                   ,&
                      & method                                                     ,&
                      & symmetry                                                   ,&
                      & coplanarity_angle                                          ,&
                      & panel_type                                                 ,&
                      & interactive                                                ,&
                      & cx,cy,cz                                                   ,&
                      & x_node,y_node,z_node                                       ,&
                      & (/l1,l2,l3/)                                               ,&
                      & (/p1,p2,p3/)                                               ,&
                      & (/n1,n2,n3/)                                               ,&
                      & speed                                                      ,&
                      & speed_x,speed_y,speed_z                                    ,&
                      & node1,node2,node3,node4                                    ,&
                      & elem1,elem2,elem3,elem4                                    ,&
                      & error                                                      )
      if (velo_err .ne. 0) then
          err = 1
          goto 999
      endif
      
      ! calculating pressure coefficients
      call pressure( interactive                                                 ,&
                  & results%cl_strip                                             ,&
                  & symmetry                                                     ,&
                  & node_num                                                     ,&
                  & panel_num                                                    ,&
                  & panel_num_no_wake                                            ,&
                  & panel_type                                                   ,&
                  & case_num_tmp                                                 ,&
                  & alfa_tmp,beta_tmp                                            ,&
                  & mac                                                          ,&
                  & wing_span                                                    ,&
                  & wing_surf                                                    ,&
                  & origin                                                       ,&
                  & speed                                                        ,&
                  & mach                                                         ,&
                  & ro                                                           ,&
                  & p_ref                                                        ,&
                  & S                                                            ,&
                  & v                                                            ,&
                  & n1,n2,n3                                                     ,&
                  & cx,cy,cz                                                     ,&
                  & x_node,y_node,z_node                                         ,&
                  & node1,node2,node3,node4                                      ,&
                  & elem1,elem2,elem3,elem4                                      ,&
                  & error                                                        ,&
                  & r_stabilizer                                                  )
      if (press_err .ne. 0) then
          err = 1
          goto 999
      endif
      
      call func_message( interactive                                              ,&
                      & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      call func_toc( interactive, start_time                                      ,&
                  & "    Time for postprocessing............................:")
      call func_new_line(interactive)
  
999 continue
  
  end subroutine calculate_postprocessing
  
  
  subroutine calculate_results(err)
      !
      implicit none
      !
      ! INTENTS OUT ==================================================================
      integer, intent(out)                                                      :: err
      ! ==============================================================================
      err = 0
      !
      ! allocate sigma to minimum size for constant doublet method
      if (method .eq. 1) allocate (sigma(1,1))
      
      ! calculating induced drag coefficient
      call induced_drag( rhs                                                      ,&
                       & sigma                                                    ,&
                       & method                                                   ,&
                       & node_num                                                 ,&
                       & panel_num                                                ,&
                       & panel_num_no_wake                                        ,&
                       & case_num                                                 ,&
                       & panel_type                                               ,&
                       & interactive                                              ,&
                       & symmetry                                                 ,&
                       & x_node,y_node,z_node                                     ,&
                       & speed                                                    ,&
                       & speed_x,speed_y,speed_z                                  ,&
                       & ro                                                       ,&
                       & wing_surf                                                ,&
                       & node1,node2,node3,node4                                  ,&
                       & elem1,elem2,elem3,elem4                                  ,& 
                       & cx,cy,cz                                                 ,&
                       & l1,l2,l3                                                 ,&
                       & p1,p2,p3                                                 ,&
                       & n1,n2,n3                                                 ,&
                       & x1,x2,x3,x4                                              ,&
                       & y1,y2,y3,y4                                              ,&
                       & d1,d2,d3,d4                                              ,&
                       & error                                                    )
      if (idrag_err .ne. 0) then
          err = 1
          goto 999
      endif
      
      deallocate( d1,d2,d3,d4                                                     ,&
                & x1,x2,x3,x4                                                     ,&
                & y1,y2,y3,y4)
      
      call func_message( interactive                                              ,&
                       & "Done!")
      call func_new_line(interactive)
      call func_new_line(interactive)
      call func_toc( interactive, start_time                                      ,&
                   & "    Time for postprocessing............................:")
      call func_new_line(interactive)
      
      ! printing coefficients
      call subr_write_coef( interactive                                           ,&
                          & case_num                                              ,&
                          & coef_x,coef_y,coef_z                                  ,&
                          & coef_l,coef_m,coef_n                                  ,&
                          & coef_drag,coef_idrag,coef_side,coef_lift              )
      
      call func_new_line(interactive)
      call func_new_line(interactive)
      
      ! WRITING APAME RESULTS FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if (res_req .eq. 1) then
      
          call func_tic(start_time)
          call func_message( interactive                                          ,&
                           & "Writing results file...")
          
          call subr_results( interactive                                          ,&
                           & version                                              ,&
                           & subversion                                           ,&
                           & build_date                                           ,&
                           & file_name                                            ,&
                           & speed                                                ,&
                           & ro                                                   ,&
                           & p_ref                                                ,&
                           & mach                                                 ,&
                           & alfa                                                 ,&
                           & beta                                                 ,&
                           & wing_span                                            ,&
                           & mac                                                  ,&
                           & wing_surf                                            ,&
                           & origin                                               ,&
                           & method                                               ,&
                           & error                                                ,&
                           & coplanarity_angle                                    ,&
                           & farfield                                             ,&
                           & collcalc                                             ,&
                           & velorder                                             ,&
                           & case_num                                             ,&
                           & node_num                                             ,&
                           & panel_num                                            ,&
                           & panel_num_no_wake                                    ,&
                           & panel_type                                           ,&
                           & x_node,y_node,z_node                                 ,&
                           & node1,node2,node3,node4                              ,&
                           & elem1,elem2,elem3,elem4                              ,&
                           & cx,cy,cz                                             ,&
                           & cp                                                   ,&
                           & v                                                    ,&
                           & vx,vy,vz                                             ,&
                           & p_dyna                                               ,&
                           & p_mano                                               ,&
                           & p_stat                                               ,&
                           & rhs                                                  ,&
                           & sigma                                                ,&
                           & S                                                    ,&
                           & FF                                                   ,&
                           & n1,n2,n3                                             ,&
                           & l1,l2,l3                                             ,&
                           & p1,p2,p3                                             ,&
                           & coef_x,coef_y,coef_z                                 ,&
                           & coef_l,coef_m,coef_n                                 ,&
                           & coef_drag,coef_idrag,coef_side,coef_lift             ,&
                           & cl_strip, cllength_strip, strip_center               ,&
                           & Fx,Fy,Fz                                             ,&
                           & Fl,Fm,Fn                                             ,&
                           & Fdrag,Fidrag,Fside,Flift                             ,&
                           & results                                              )
          
          call func_message( interactive                                          ,&
                           & "Done!")
          call func_new_line(interactive)
          call func_new_line(interactive)
          call func_toc( interactive, start_time                                  ,&
                       & "    Time for writing results file......................:")
          call func_new_line(interactive)
      endif
      
      if (vtk_req .eq. 1) then
      
          ! allocate sigma to minimum size for constant doublet method
          if ((method .eq. 1) .and. (allocated(sigma) .eq. .false.)) allocate (sigma(1,1))
      
          call func_tic(start_time)
          call func_message( interactive                                          ,&
                           & "Writing vtk file...")
      
          call subr_vtk_results( interactive                                      ,&
                           & version                                              ,&
                           & subversion                                           ,&
                           & build_date                                           ,&
                           & file_name                                            ,&
                           & speed                                                ,&
                           & ro                                                   ,&
                           & p_ref                                                ,&
                           & mach                                                 ,&
                           & alfa                                                 ,&
                           & beta                                                 ,&
                           & wing_span                                            ,&
                           & mac                                                  ,&
                           & wing_surf                                            ,&
                           & origin                                               ,&
                           & method                                               ,&
                           & error                                                ,&
                           & coplanarity_angle                                    ,&
                           & farfield                                             ,&
                           & collcalc                                             ,&
                           & velorder                                             ,&
                           & case_num                                             ,&
                           & node_num                                             ,&
                           & panel_num                                            ,&
                           & panel_num_no_wake                                    ,&
                           & panel_type                                           ,&
                           & x_node,y_node,z_node                                 ,&
                           & node1,node2,node3,node4                              ,&
                           & elem1,elem2,elem3,elem4                              ,&
                           & cx,cy,cz                                             ,&
                           & cp                                                   ,&
                           & v                                                    ,&
                           & vx,vy,vz                                             ,&
                           & p_dyna                                               ,&
                           & p_mano                                               ,&
                           & p_stat                                               ,&
                           & rhs                                                  ,&
                           & sigma                                                ,&
                           & S                                                    ,&
                           & FF                                                   ,&
                           & n1,n2,n3                                             ,&
                           & l1,l2,l3                                             ,&
                           & p1,p2,p3                                             ,&
                           & coef_x,coef_y,coef_z                                 ,&
                           & coef_l,coef_m,coef_n                                 ,&
                           & coef_drag,coef_idrag,coef_side,coef_lift             ,&
                           & cl_strip, cllength_strip, elem_strip_id              ,&
                           & Fx,Fy,Fz                                             ,&
                           & Fl,Fm,Fn                                             ,&
                           & Fdrag,Fidrag,Fside,Flift                             ,&
                           & results                                              )
      
          call func_message( interactive                                          ,&
                           & "Done!")
          call func_new_line(interactive)
          call func_new_line(interactive)
          call func_toc( interactive, start_time                                  ,&
                       & "    Time for writing vtk file......................:")
          call func_new_line(interactive)
      endif
      
      write(*,*) 'coef_idrag:   ', coef_idrag
      write(*,*) 'coef_lift:    ', coef_lift 
      
      call func_message( interactive                                              ,&
                       & "Job complete!")
      call func_new_line(interactive)
  
999 continue
  
  end subroutine calculate_results


  subroutine deallocate_arrays()
      !
      implicit none
      !
      ! deallocate arrays
      if (allocated(speed_x)) deallocate(speed_x)
      if (allocated(speed_y)) deallocate(speed_y)
      if (allocated(speed_z)) deallocate(speed_z)
      if (allocated(sigma)) deallocate(sigma)
      if (allocated(rhs)) deallocate(rhs)
      if (allocated(v)) deallocate(v)
      if (allocated(vx)) deallocate(vx)
      if (allocated(vy)) deallocate(vy)
      if (allocated(vz)) deallocate(vz)
      if (allocated(cp)) deallocate(cp)
      if (allocated(p_dyna)) deallocate(p_dyna)
      if (allocated(p_mano)) deallocate(p_mano)
      if (allocated(p_stat)) deallocate(p_stat)
      if (allocated(Ftot)) deallocate(Ftot)
      if (allocated(Fx)) deallocate(Fx)
      if (allocated(Fy)) deallocate(Fy)
      if (allocated(Fz)) deallocate(Fz)
      if (allocated(Fl)) deallocate(Fl)
      if (allocated(Fm)) deallocate(Fm)
      if (allocated(Fn)) deallocate(Fn)
      if (allocated(Fdrag)) deallocate(Fdrag)
      if (allocated(Fside)) deallocate(Fside)
      if (allocated(Flift)) deallocate(Flift)
      if (allocated(coef_x)) deallocate(coef_x)
      if (allocated(coef_y)) deallocate(coef_y)
      if (allocated(coef_z)) deallocate(coef_z)
      if (allocated(coef_l)) deallocate(coef_l)
      if (allocated(coef_m)) deallocate(coef_m)
      if (allocated(coef_n)) deallocate(coef_n)
      if (allocated(coef_drag)) deallocate(coef_drag)
      if (allocated(coef_side)) deallocate(coef_side)
      if (allocated(coef_lift)) deallocate(coef_lift)
      if (allocated(cl_mom_corrected)) deallocate(cl_mom_corrected)
      if (allocated(Fidrag)) deallocate(Fidrag)
      if (allocated(coef_idrag)) deallocate(coef_idrag)
      if (allocated(cl_strip)) deallocate(cl_strip)
      if (allocated(cllength_strip)) deallocate(cllength_strip)
      if (allocated(strip_center)) deallocate(strip_center)
      if (allocated(elem_strip_id)) deallocate(elem_strip_id)
      !
  end subroutine deallocate_arrays
  
  subroutine deallocate_grid()
      !
      implicit none
      !
      ! deallocate array
      if (allocated(a)) deallocate(a)
      if (allocated(a_no_wake)) deallocate(a_no_wake)
      if (allocated(b)) deallocate(b)
      if (allocated(FF)) deallocate(FF)
      if (allocated(S)) deallocate(S)
      if (allocated(n1)) deallocate(n1)
      if (allocated(n2)) deallocate(n2)
      if (allocated(n3)) deallocate(n3)
      if (allocated(l1)) deallocate(l1)
      if (allocated(l2)) deallocate(l2)
      if (allocated(l3)) deallocate(l3)
      if (allocated(p1)) deallocate(p1)
      if (allocated(p2)) deallocate(p2)
      if (allocated(p3)) deallocate(p3)
      if (allocated(x1)) deallocate(x1)
      if (allocated(x2)) deallocate(x2)
      if (allocated(x3)) deallocate(x3)
      if (allocated(x4)) deallocate(x4)
      if (allocated(y1)) deallocate(y1)
      if (allocated(y2)) deallocate(y2)
      if (allocated(y3)) deallocate(y3)
      if (allocated(y4)) deallocate(y4)
      if (allocated(d1)) deallocate(d1)
      if (allocated(d2)) deallocate(d2)
      if (allocated(d3)) deallocate(d3)
      if (allocated(d4)) deallocate(d4)
      if (allocated(cx)) deallocate(cx)
      if (allocated(cy)) deallocate(cy)
      if (allocated(cz)) deallocate(cz)
      if (allocated(colx)) deallocate(colx)
      if (allocated(coly)) deallocate(coly)
      if (allocated(colz)) deallocate(colz)
      !
  end subroutine deallocate_grid

end subroutine apame_main
