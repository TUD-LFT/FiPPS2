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
subroutine testcase1 (testcase)
			
! =================================================================================================
!
!	Header:		subroutine defining testcase1
!
!	Content:	1 tria3-element in global coordinates
!			1 force in y-direction on node 3
!			1 x,y,z,rotx,roty-bc on nodes 1,2
!			
!
!	Input:		
!
!	Output:		
!
!	Calls:		
!
!	Called by:	
!
!	Author:		Martin Rädel			08.12.2009
! 			TU Dresden, Diplomarbeit
!
! =================================================================================================
!
use globale_variablen
!
! =================================================================================================
!
implicit none
!
! =================================================================================================
!
! Input
!
integer, intent(in)	:: testcase
!
! Output
!
!
! Input + Output
!




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

!
! =================================================================================================

! 1 tria3 element in global system, in-plane tension on node 3
         
is_node    = .true.
is_tria3   = .true.
is_load    = .true.
is_force   = .true.
is_mat1    = .true.
is_pshell  = .true.
is_spcadd  = .true.
is_spc1    = .true.
is_subcase = .true.
	 
! nodes
	 
allocate(nodes(1:3))
	 
!	 nodes(1) = (/ 1,0,0.D0 ,0.D0 ,0.D0  /)
!	 nodes(2) = (/ 2,0,1.D0 ,0.D0 ,0.D0  /)
!	 nodes(3) = (/ 3,0,0.5D0,1.D0 ,0.D0  /)
         
nodes(1)%nid         = 1
nodes(1)%cid         = 0
nodes(1)%coords(1:3) = (/ 0.D0, 0.D0, 0.D0 /)

nodes(2)%nid         = 2
nodes(2)%cid         = 0
nodes(2)%coords(1:3) = (/ 1.D0, 0.D0, 0.D0 /)

nodes(3)%nid         = 3
nodes(3)%cid         = 0
nodes(3)%coords(1:3) = (/ 0.5D0,1.D0, 0.D0 /)

! tria3 elements

allocate(tria3s(1))
	 
!	 tria3s(1) = (/ 1,1,1,2,3 /)

tria3s(1)%eid       = 1
tria3s(1)%pid       = 1
tria3s(1)%nids(1:3) = (/ 1,2,3 /)
tria3s(1)%theta     = 0.D0
tria3s(1)%offset    = 0.D0
	 
! loads
	 
allocate(loads(1))
	 
!	 loads(1) = (/ 1,1.D0, 1.D0 ,1 /)
	 
loads(1)%lcid  = 1
loads(1)%sfac  = 1.D0
loads(1)%sfaci = 1.D0
loads(1)%lidi  = 1
	 
! forces
	 
allocate(forces(1))
	 
!	 forces(1) = (/ 1, 3, 0, 10.D0 , 0.D0 , 1.D0 , 0.D0 /)
	 
forces(1)%lid = 1
forces(1)%nid = 3
forces(1)%cid = 0
forces(1)%fac = 10.D0
forces(1)%ni  = (/ 0.D0, 1.D0, 0.D0 /)
	 
! pshells

allocate(pshells(1))
	 
!	 pshells(1) = (/ 1, 1, 1.D0, 1, 1.D0, 1, 0.833333333D0, 0.D0, -0.5D0, 0.5D0, 1 /)
	 
pshells(1)%pid  = 1
pshells(1)%mid1 = 1
pshells(1)%mt   = 1.D0
pshells(1)%mid2 = 1
pshells(1)%bmr  = 1.D0
pshells(1)%mid3 = 1
pshells(1)%tst  = 0.8333333333D0
pshells(1)%nsm  = 0.D0
pshells(1)%z1 	 = -0.5D0
pshells(1)%z2   =  0.5D0
pshells(1)%mid4 = 1

! mat1s

allocate(mat1s(1))

!	 mat1s(1) = (/ 1, 210000.D0, 78947.37D0, 0.33D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
	 
mat1s(1)%mid  = 1
mat1s(1)%ym   = 210000.D0
mat1s(1)%sm   = 78947.37D0
mat1s(1)%nu   = 0.33D0
mat1s(1)%rho  = 0.D0
mat1s(1)%ath  = 0.D0
mat1s(1)%tref = 0.D0
mat1s(1)%ge   = 0.D0

! spcadds

allocate(spcadds(1))

!	 spcadds(1) = (/ 1, 1 /)
	 
spcadds(1)%scid = 1
spcadds(1)%sid  = 1

! spc1s

allocate(spc1s(1))

!	 spc1s(1) = (/ 1, 12345, 1, .true., 2 /)

spc1s(1)%sid = 1
spc1s(1)%dof = 12345
spc1s(1)%n1  = 1
spc1s(1)%thru= .true.
spc1s(1)%nn  = 2

! subcases

allocate(subcases(1))

!	 subcases(1) = (/ 1, 1, 1 /)
	 
subcases(1)%scid     = 1
subcases(1)%spcaddid = 1
subcases(1)%loadid   = 1

! solution 1-static, 2-static+buckling

sol=1

end subroutine testcase1
