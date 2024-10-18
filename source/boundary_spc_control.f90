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
!> control subroutine for inserting homogeneous boundary conditions
!
!> @details
!
!> @author Andreas Hauffe, TU Dresden, wissenschaftlicher Mitarbeiter 10.10.2012
!
!> $Id: boundary_spc_control.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine boundary_spc_control (fesim,fixedU,scloop)

use konstanten
use fesimulation_typen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
type(fe_simulation)                     :: fesim
integer,intent(in)                      :: scloop
!
! Output
!
logical, dimension(fesim%num_dof),intent(out) :: fixedU
!
! inner
!
integer                                 :: count_dof
integer                                 :: gg,ii,jj,kk
integer                                 :: err_code=0
integer, dimension(:), allocatable      :: mod_dof
logical, dimension(:), allocatable      :: help

! Loop over homogeneous bc

if (fesim%is_spcadd == .true.) then
   do gg=1,size(fesim%randbedingungen%spcadds,1)
      if (fesim%lasten%subcases(scloop)%spcaddid == fesim%randbedingungen%spcadds(gg)%scid) then
         
         ! Single point constraints
         
         if (fesim%is_spc1 == .true.) then
            do ii=1,size(fesim%randbedingungen%spc1s,1)
               if (fesim%randbedingungen%spcadds(gg)%sid == fesim%randbedingungen%spc1s(ii)%sid) then
                  
                  ! calculate vector of dofs size
                  
                  if (fesim%randbedingungen%spc1s(ii)%dof > 10d4) then

                     count_dof = 6
                     
                  else if (fesim%randbedingungen%spc1s(ii)%dof > 10d3) then
                        
                     count_dof = 5

                  else if (fesim%randbedingungen%spc1s(ii)%dof > 10d2) then
                     
                     count_dof = 4
                     
                  else if (fesim%randbedingungen%spc1s(ii)%dof > 10d1) then
                     
                     count_dof = 3
                     
                  else if (fesim%randbedingungen%spc1s(ii)%dof > 10d0) then
                     
                     count_dof = 2
                     
                  else if (fesim%randbedingungen%spc1s(ii)%dof > 0) then
                     
                     count_dof = 1

                  else
   
                    write(*,*) 'Error reading spc1s value #dof# on line', ii
                    err_code=1
                    goto 9999
                     
                  end if

                  ! allocate memory
                  
                  allocate(mod_dof(count_dof))
                  
                  ! get dofs

                  call spc1_mod_dof(count_dof,fesim%randbedingungen%spc1s(ii)%dof,mod_dof)
                  
                  ! loop over nodes
                  
                  if (fesim%randbedingungen%spc1s(ii)%thru == .false.) then             ! bc on 1 node only

                    do kk = 1, count_dof
                      fixedU((fesim%randbedingungen%spc1s(ii)%n1-1)*ndof+mod_dof(kk)) = .TRUE.
                    end do
                  
                  else if (fesim%randbedingungen%spc1s(ii)%thru == .true.) then ! bc on multiple nodes
                  
                    do jj=fesim%randbedingungen%spc1s(ii)%n1,fesim%randbedingungen%spc1s(ii)%nn

                      do kk = 1, count_dof
                        fixedU((jj-1)*ndof+mod_dof(kk)) = .TRUE.
                      end do
                        
                    end do
                  
                  else
                     
                    write(*,*) 'Error reading spc1s value #thru# on line', ii
                    err_code=1
                    goto 9999
                    
                  end if

                  deallocate(mod_dof)
                  
               end if
            end do
         end if

    ! Enforced motion constraints
 
    if (fesim%is_spcd == .true.) then
      do ii=1,size(fesim%randbedingungen%spcds,1)
        if (fesim%randbedingungen%spcadds(gg)%sid == fesim%randbedingungen%spcds(ii)%sid) then
  
       ! set dof as constraint
  
       ! loop over nodes

       ! do kk=1,size(nodes,1)

         ! if (fesim%randbedingungen%spcds(ii)%nid == nodes(kk)%nid) then

           ! check for constraint

            if (fixedU((fesim%randbedingungen%spcds(ii)%nid-1)*ndof+fesim%randbedingungen%spcds(ii)%dof) == .TRUE.) then

              write(*,*) 'FiPPS WARNING:'
              write(*,*) 'you try to set an enforced motion constraint on a dof'
              write(*,*) 'that already is locked via a spc1 constraint'
              write(*,*) 'enforced motion constraint will have no effect'

            else

              fixedU((fesim%randbedingungen%spcds(ii)%nid-1)*ndof+fesim%randbedingungen%spcds(ii)%dof) = .TRUE.

            end if

           ! end if

         ! end do

         end if
       end do
     end if
    end if
  end do
end if

! If solid elements are used, the rotational degrees of freedom
! of all nodes, not part of any other beam or shell element have
! to be fixed. The solid element only takes translational degrees
! of freedom into account.
if (fesim%is_lsolid20 == .true.) then

  allocate(help(fesim%num_nodes))
  
  ! Nodes not part of any beam or shell element get false
  help(:) = .false.

  ! Loop over all nodes of beam elements
  if (fesim%is_beam2 == .true.) then
    do ii=1,size(fesim%elemente%beam2s,1)
      do jj=1,2
        help(fesim%elemente%beam2s(ii)%nids(jj)) = .true.
      end do
    end do
  end if

  ! Loop over all nodes of shell elements
  if (fesim%is_quad8 == .true.) then
    do ii=1,size(fesim%elemente%quad8s,1)
      do jj=1,8
        help(fesim%elemente%quad8s(ii)%nids(jj)) = .true.
      end do
    end do
  end if

  ! Fix rotational degrees of freedom of all remaining nodes
  do ii=1,fesim%num_nodes
    if (help(ii) .eqv. .false.) then
      fixedU((ii-1)*ndof+4) = .TRUE.
      fixedU((ii-1)*ndof+5) = .TRUE.
      fixedU((ii-1)*ndof+6) = .TRUE.
    end if
  end do
  
  deallocate(help)

  write(*,*) 'Rotationen fuer Lsolid20-Elemente festgesetzt!'
end if

!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'boundary_spc_control'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

return

end subroutine boundary_spc_control
