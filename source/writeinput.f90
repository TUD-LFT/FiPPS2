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
!> $Id: writeinput.f90 484 2024-10-18 14:28:29Z s1080304 $
!> $Author: s1080304 $
!> $Revision: 484 $
!> $Date: 2024-10-18 16:28:29 +0200 (Fr, 18. Okt 2024) $
!
! =================================================================================================
subroutine writeInput()

  use globale_variablen
  use netz_variablen
!
! =================================================================================================
!
  implicit none
!
! =================================================================================================
!
! Include
!

!
! =================================================================================================
!
! Data types
!  
  integer :: out = 20
  integer :: ii, jj
  
  integer       :: err_code=0
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
  open(unit=out, file='input.log', status='UNKNOWN') 
  
  write(out,*) 'Nodes'
  if (is_node == .true.) then
    do ii = 1, size(nodes,1)
      write(out,'(3E24.16)') (nodes(ii)%coords(jj), jj=1,3)
    end do
  end if
  write(out,*) 'Loads'
  if (is_load == .true.) then
    do ii = 1, size(loads,1)
      write(out,'(I10,E23.16,I10)') loads(ii)%lcid, loads(ii)%sfaci, loads(ii)%lidi
    end do
  end if
  
  write(out,*) 'SPCAdd'
  if (is_spcadd == .true.) then
    do ii = 1, size(spcadds,1)
      write (out,'(2I10)') spcadds(ii)%scid, spcadds(ii)%sid
    end do
  end if
  
  write(out,*) 'Beam2'
  if (is_beam2 == .true.) then
    do ii = 1, size(beam2s,1)
      write (out,'(3I10,3E23.16,I10)') beam2s(ii)%pid,(beam2s(ii)%nids(jj),jj=1,2), (beam2s(ii)%xi(jj), jj=1,3), beam2s(ii)%n0
    end do
  end if
  
  write(out,*) 'Tria3'
  if (is_tria3 == .true.) then
    do ii = 1, size(tria3s,1)
      write (out,'(4I10,E23.16)') tria3s(ii)%pid,(tria3s(ii)%nids(jj),jj=1,3), tria3s(ii)%theta
    end do
  end if
  
  write(out,*) 'Quad4'
  if (is_quad4 == .true.) then
    do ii = 1, size(quad4s,1)
      write (out,'(5I10,E23.16)') quad4s(ii)%pid,(quad4s(ii)%nids(jj),jj=1,4), quad4s(ii)%theta
    end do
  end if
  
  write(out,*) 'Quad8'
  if (is_quad8 == .true.) then
    do ii = 1, size(quad8s,1)
      write (out,'(9I10,E23.16)') quad8s(ii)%pid,(quad8s(ii)%nids(jj),jj=1,8), quad8s(ii)%theta
    end do
  end if
  
  write(out,*) 'LSolid20'
  if (is_lsolid20 == .true.) then
    do ii = 1, size(lsolid20s,1)
      write (out,'(21I10)') lsolid20s(ii)%pid,(lsolid20s(ii)%nids(jj),jj=1,20)
    end do
  end if
  
  write(out,*) 'Force'
  if (is_force == .true.) then
    do ii = 1, size(forces,1)
      write (out,'(2I10,4E23.16)') forces(ii)%lid, forces(ii)%nid, forces(ii)%fac, (forces(ii)%ni(jj), jj=1,3)
    end do
  end if
  
  write(out,*) 'Moment'
  if (is_moment == .true.) then
    do ii = 1, size(moments,1)
      write (out,'(2I10,4E23.16)') moments(ii)%lid,moments(ii)%nid,moments(ii)%fac,(moments(ii)%ni(jj), jj=1,3)
    end do
  end if
  
  write(out,*) 'P3Load'
  if (is_p3load == .true.) then
    do ii = 1, size(p3loads,1)
      write (out,'(3I10,3E23.16,L,I10)') p3loads(ii)%lid,p3loads(ii)%eid1,p3loads(ii)%cid,(p3loads(ii)%pi(jj), jj=1,3), &
                                                                                & p3loads(ii)%thru, p3loads(ii)%eid2
    end do
  end if

  write(out,*) 'P2Load'
  if (is_p2load == .true.) then
    do ii = 1, size(p2loads,1)
      write (out,'(3I10,2E23.16,L1,I10)') p2loads(ii)%lid,p2loads(ii)%eid1,p2loads(ii)%dir, &
                                        & (p2loads(ii)%pi(jj), jj=1,2),p2loads(ii)%thru,p2loads(ii)%eid2
    end do
  end if

  write(out,*) 'P8Load'
  if (is_p8load == .true.) then
    do ii = 1, size(p8loads,1)
      write (out,'(3I10,4E23.16,L1,I10)') p8loads(ii)%lid,p8loads(ii)%eid1,p8loads(ii)%cid,(p8loads(ii)%pi(jj), jj=1,4), &
                                        & p8loads(ii)%thru,p8loads(ii)%eid2
    end do
  end if

  write(out,*) 'P20Load'
  if (is_p20load == .true.) then
    do ii = 1, size(p20loads,1)
      write (out,'(2I10,E23.16,I10,L1,I10)') p20loads(ii)%lid,p20loads(ii)%eid1,p20loads(ii)%p,p20loads(ii)%surf, &
                                           & p20loads(ii)%thru,p20loads(ii)%eid2
    end do
  end if

  write(out,*) 'Beam2Temp'
  if (is_beam2temp == .true.) then
    do ii = 1, size(beam2temps,1)
      write (out,'(2I10,E23.16)') beam2temps(ii)%lid,beam2temps(ii)%eid,beam2temps(ii)%temp
    end do
  end if

  write(out,*) 'Quad8Temp'
  if (is_quad8temp == .true.) then
    do ii = 1, size(quad8temps,1)
      write (out,'(2I10,E23.16)') quad8temps(ii)%lid,quad8temps(ii)%eid,quad8temps(ii)%temp
    end do
  end if

  write(out,*) 'Lsolid20Temp'
  if (is_lsolid20temp == .true.) then
    do ii = 1, size(lsolid20temps,1)
      write (out,'(2I10,E23.16)') lsolid20temps(ii)%lid,lsolid20temps(ii)%eid,lsolid20temps(ii)%temp
    end do
  end if

  write(out,*) 'Aeroload2d'
  if (is_aeroload2d == .true.) then
    do ii = 1, size(aeroload2ds,1)
      write (out,'(2I10)') aeroload2ds(ii)%lid,aeroload2ds(ii)%mthd
    end do
  end if

  write(out,*) 'Aeroload3d'
  if (is_aeroload3d == .true.) then
    do ii = 1, size(aeroload3ds,1)
      write (out,'(2I10)') aeroload3ds(ii)%lid,aeroload3ds(ii)%mthd
    end do
  end if
  
  write(out,*) 'Mat1'
  if (is_mat1 == .true.) then
    do ii = 1, size(mat1s,1)
      write (out,'(I10,7E23.16)') mat1s(ii)%mid, mat1s(ii)%ym, mat1s(ii)%sm, mat1s(ii)%nu, mat1s(ii)%rho, mat1s(ii)%ath,&
                                & mat1s(ii)%tref, mat1s(ii)%ge
    end do
  end if
  
  write(out,*) 'Mat8'
  if (is_mat8 == .true.) then
    do ii = 1, size(mat8s,1)
      write (out,'(I10,11E23.16)') mat8s(ii)%mid,mat8s(ii)%ym11,mat8s(ii)%ym22,mat8s(ii)%nu12,mat8s(ii)%sm12,mat8s(ii)%sm13,&
                                 & mat8s(ii)%sm23,mat8s(ii)%rho,mat8s(ii)%ath11,mat8s(ii)%ath22,mat8s(ii)%tref,mat8s(ii)%ge
    end do
  end if
  
  write(out,*) 'Mat20'
  if (is_mat20 == .true.) then
    do ii = 1, size(mat20s,1)
      write (out,'(I10,13E23.16)') mat20s(ii)%mid,mat20s(ii)%ym11,mat20s(ii)%ym22,mat20s(ii)%ym33,mat20s(ii)%nu12,mat20s(ii)%nu13, &
                                    & mat20s(ii)%nu23,mat20s(ii)%sm12,mat20s(ii)%sm13,mat20s(ii)%sm23,mat20s(ii)%ath11,mat20s(ii)%ath22,mat20s(ii)%ath33
    end do
  end if
  
  write(out,*) 'PShell'
  if (is_pshell == .true.) then
    do ii = 1, size(pshells,1)
      write (out,'(2I10,E23.16,I10,E23.16,I10,4E23.16,I10)') pshells(ii)%pid, pshells(ii)%mid1, pshells(ii)%mt, &
                                        & pshells(ii)%mid2, pshells(ii)%bmr, pshells(ii)%mid3, pshells(ii)%tst, pshells(ii)%nsm, pshells(ii)%z1, pshells(ii)%z2, &
                                        & pshells(ii)%mid4
    end do
  end if
  
  write(out,*) 'PComps'
  if (is_pcomp == .true.) then
    do ii = 1, size(pcomps,1)
      write (out,'(2I10,E23.16,I10,2E23.16,A2)') pcomps(ii)%pid,pcomps(ii)%lamid,pcomps(ii)%offset,pcomps(ii)%lay,&
                                                                                & pcomps(ii)%nsm,pcomps(ii)%sb,pcomps(ii)%ft
    end do
  end if
  
  write(out,*) 'Plsolid'
  if (is_plsolid == .true.) then
    do ii = 1, size(plsolids,1)
      write (out,'(3I10,L1,I10)') plsolids(ii)%pid,plsolids(ii)%lamid,plsolids(ii)%cid,plsolids(ii)%globOut,plsolids(ii)%resLay
    end do
  end if
  
  write(out,*) 'SPC1'
  if (is_spc1 == .true.) then
    do ii = 1, size(spc1s,1)
      write (out,'(3I10,L1,I10)') spc1s(ii)%sid, spc1s(ii)%dof, spc1s(ii)%n1, spc1s(ii)%thru, spc1s(ii)%nn
    end do
  end if
  
  write(out,*) 'Coord'
  if (is_coord == .true.) then
    do ii = 1, size(coords,1)
      write(out,'(I10,9E23.16)') coords(ii)%cid, (coords(ii)%xAxisVec(jj), jj=1,3), &
                                                 (coords(ii)%yAxisVec(jj), jj=1,3), &
                                                 (coords(ii)%zAxisVec(jj), jj=1,3)
    end do
  end if
  
  write(out,*) 'Subcase'
  if (is_subcase == .true.) then
    do ii = 1, size(subcases,1)
      write (out,'(3I10)') subcases(ii)%scid, subcases(ii)%spcaddid, subcases(ii)%loadid
    end do
  end if
  
  write(out,*) 'Lam8'
  if (is_lam8 == .true.) then
    do ii = 1, size(lam8s,1)
      write (out,'(3I10,2E23.16,A3)') lam8s(ii)%lamid,lam8s(ii)%plyid,lam8s(ii)%mat8id,lam8s(ii)%th,&
                                                                        & lam8s(ii)%angle,'rad'
    end do
  end if
  
  write(out,*) 'Lam20'
  if (is_lam20 == .true.) then
    do ii = 1, size(lam20s,1)
      write (out,'(3I10,2E23.16,A3,I10)') lam20s(ii)%lamid,lam20s(ii)%plyid,lam20s(ii)%mat20id(msid),lam20s(ii)%th,&
                                        & lam20s(ii)%angle,'rad',lam20s(ii)%nop
    end do
  end if
  
  write(out,*) 'Coupling'
  if (is_coupling == .true.) then
    do ii = 1, size(couplings,1)
      write (out,'(3I10)') couplings(ii)%cpsid,couplings(ii)%dof,couplings(ii)%nid
    end do
  end if
  
  close(out)
  
!
! =================================================================================================
!
! Error handling
!
9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)                      'writeInput'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if
  
end subroutine writeInput
