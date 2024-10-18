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
subroutine testcase12 (testcase)
			
! =================================================================================================
!
!	Header:		subroutine defining testcase12
!
!	Content:	1 flat quad8-element in global coordinates
!			1 force in y-direction on node 3
!			1 x,y,z,rotx,roty-bc on nodes 1,2,3
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
integer 		:: io_error, read_error			! Error-handling parameters
integer 		:: rowcount
integer			:: err_code=0

integer			:: num
integer			:: ii,jj

 character(20)		:: fname

!	 __	 _______________
!	 |	|		|
!	 b	|		|
!	 r	|		|
!	 e 	|		|
!  	 i	|		|
!	 t	|		|
!  yA	 e	|		|
!   |	 |_	|_______________|
!
!		|<---laenge---->|
!
!		->x
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
is_quad8   = .true.
is_load    = .true.
is_force   = .true.
is_mat1    = .true.
is_pshell  = .true.
is_mat8    = .true.
is_pcomp   = .true.
is_lam8    = .true.
is_spcadd  = .true.
is_spc1    = .true.
is_subcase = .true.


! NODES

if (is_node == .false.) then
	write(*,*) 'No nodes defined on ', fname
	err_code=1
	goto 9999
else
	fname='nodes'
	open (unit=20,file=fname,status='old',action='read', iostat=io_error)
		! unit=x, x=Dateinummer (Ganzzahl,>10,da oft Nummern <10 fix zugeordnet, z.B. der Standardein-, ausgabe)
		! file=x, x=Dateiname
		! status, old-Datei existiert bereits. default: 'unknown'
		! action, read-lesen, write-schreiben, readwrite-lesen+schreiben
	if (io_error == 0) then
		! read 1st actual value
		read (20,*,iostat=read_error) rowcount	! Zeilenanzahl aus erster Zeile
		! error handling
		if (read_error /= 0) then
			write(*,*) 'Error reading file ', fname
			write(*,*) 'number of rows could not be determined from first row'
			err_code=1
			goto 9999
		end if
		! allocate memory for type and number of rows
		allocate(nodes(1:rowcount))
		! set fixed values
		do ii=1,rowcount
		   nodes(ii)%nid = ii
		   nodes(ii)%cid = 0
		end do
		! read values
		do ii=1,rowcount
!			read (20,'(3E24.16)',iostat=read_error) (nodes(ii)%coords(jj), jj=1,3)
			read (20,'(I10,X,3E24.16)',iostat=read_error) nodes(ii)%nid,(nodes(ii)%coords(jj), jj=1,3)
			if (read_error /= 0) then
				write(*,*) 'Error reading file ', fname
				err_code=1
				goto 9999
			end if
		end do
		! Close file
		close (20, iostat=io_error)
		! Error handling
		if (io_error /= 0) then
			write(*,*) 'Error closing file ', fname
			err_code=1
			goto 9999
		end if
	else
		write (*,*) 'Error opening file ', fname
		err_code=1
		goto 9999
	end if
end if

! ELEMENTS

if (is_quad8 /= .false.) then
	fname='elements'
	open (unit=37,file=fname,status='old',action='read', iostat=io_error)
	if (io_error == 0) then
		read (37,*,iostat=read_error) rowcount
		if (read_error /= 0) then
			write(*,*) 'Error reading file ', fname
			err_code=1
			goto 9999
		end if
		allocate(quad8s(1:rowcount))
		! set fixed values
		do ii=1,rowcount
		   quad8s(ii)%eid    = ii
		   quad8s(ii)%pid    = 2
		   quad8s(ii)%theta  = 0.D0
		   quad8s(ii)%atype  = 'deg'
		   quad8s(ii)%offset = 0.D0
		end do
		! read values
		do ii=1,rowcount
			read (37,'(8I8)',iostat=read_error) (quad8s(ii)%nids(jj),jj=1,8)
			if (read_error /= 0) then
				write(*,*) 'Error reading file ', fname
				err_code=1
				goto 9999
			end if
		end do
		close (37, iostat=io_error)
		if (io_error /= 0) then
			write(*,*) 'Error closing file ', fname
			err_code=1
			goto 9999
		end if
	else
		write (*,*) 'Error opening file ', fname
		err_code=1
		goto 9999
	end if
end if

! FORCES

if (is_force /= .false.) then
	fname='forces'
	open (unit=25,file=fname,status='old',action='read', iostat=io_error)
	if (io_error == 0) then
		read (25,*,iostat=read_error) rowcount
		if (read_error /= 0) then
			write(*,*) 'Error reading file ', fname
			err_code=1
			goto 9999
		end if
		allocate(forces(1:rowcount))
		! set fixed values
		do ii=1,rowcount
		   forces(ii)%lid = ii
		   forces(ii)%cid = 0
		   forces(ii)%fac = -1.D0
		end do
		! read values		
		do ii=1,rowcount
			read (25,'(I8,3E24.16)',iostat=read_error) forces(ii)%nid, (forces(ii)%ni(jj), jj=1,3)
			if (read_error /= 0) then
				write(*,*) 'Error reading file ', fname
				err_code=1
				goto 9999
			end if
		end do
		close (25, iostat=io_error)
		if (io_error /= 0) then
			write(*,*) 'Error closing file ', fname
			err_code=1
		goto 9999
		end if
	else
		write (*,*) 'Error opening file ', fname
		err_code=1
		goto 9999
	end if
end if

! DISPLACEMENTS

if (is_spc1 /= .false.) then
	fname='displacements'
	open (unit=33,file=fname,status='old',action='read', iostat=io_error)
	if (io_error == 0) then
		read (33,*,iostat=read_error) rowcount
		if (read_error /= 0) then
			write(*,*) 'Error reading file ', fname
			err_code=1
			goto 9999
		end if
		allocate(spc1s(1:rowcount))
		! set fixed values
		do ii=1,rowcount
		   spc1s(ii)%sid  = ii
		   spc1s(ii)%thru = .FALSE.
		   spc1s(ii)%nn   = 0
		end do
		! read values
		do ii=1,rowcount
			read (33,'(2I10)',iostat=read_error) spc1s(ii)%n1, spc1s(ii)%dof
			if (read_error /= 0) then
				write(*,*) 'Error reading file ', fname
				err_code=1
				goto 9999
			end if
		end do
		close (33, iostat=io_error)
		if (io_error /= 0) then
			write(*,*) 'Error closing file ', fname
			err_code=1
			goto 9999
		end if
	else
		write (*,*) 'Error opening file ', fname
		err_code=1
		goto 9999
	end if
end if
	 
! loads
	 
allocate(loads(size(forces,1)))

do ii=1,size(forces,1)
	 
loads(ii)%lcid  = 1
loads(ii)%sfac  = 1.D0
loads(ii)%sfaci = 1.D0
loads(ii)%lidi  = ii

end do
	 
! pshells

allocate(pshells(1))
        
pshells(1)%pid  = 1
pshells(1)%mid1 = 1
pshells(1)%mt	= 2.D0
pshells(1)%mid2 = 1
pshells(1)%bmr  = 1.D0
pshells(1)%mid3 = 1
pshells(1)%tst  = 0.8333333333D0
pshells(1)%nsm  = 0.D0
pshells(1)%z1	 = -0.5D0
pshells(1)%z2	=  0.5D0
pshells(1)%mid4 = 1

! mat1s

allocate(mat1s(1))
        
!mat1s(1)%mid  = 1
!mat1s(1)%ym   = 210000.D0
!mat1s(1)%sm   = 78947.37D0
!mat1s(1)%nu   = 0.33D0
!mat1s(1)%rho  = 0.D0
!mat1s(1)%ath  = 0.D0
!mat1s(1)%tref = 0.D0
!mat1s(1)%ge   = 0.D0

mat1s(1)%mid  = 1
mat1s(1)%ym   = 70000.D0
mat1s(1)%sm   = 26923.07D0
mat1s(1)%nu   = 0.3D0
mat1s(1)%rho  = 0.D0
mat1s(1)%ath  = 0.D0
mat1s(1)%tref = 0.D0
mat1s(1)%ge   = 0.D0

! pcomps

allocate(pcomps(1))

pcomps(1)%pid	 = 2
pcomps(1)%lamid  = 1
pcomps(1)%offset = 0.D0
pcomps(1)%lay	 = 8
pcomps(1)%nsm	 = 0
pcomps(1)%sb	 = 0.D0

! lam8s

allocate(lam8s(8))

! 1. Schicht

lam8s(1)%lamid  = 1
lam8s(1)%plyid  = 1
lam8s(1)%mat8id = 2
lam8s(1)%th	= 0.25D0
lam8s(1)%angle  = 0.D0
lam8s(1)%atype  = 'deg'

! 2. Schicht

lam8s(2)%lamid  = 1
lam8s(2)%plyid  = 2
lam8s(2)%mat8id = 2
lam8s(2)%th	= 0.25D0
lam8s(2)%angle  = 90.D0
lam8s(2)%atype  = 'deg'

! 3. Schicht

lam8s(3)%lamid  = 1
lam8s(3)%plyid  = 3
lam8s(3)%mat8id = 2
lam8s(3)%th	= 0.25D0
lam8s(3)%angle  = 0.D0
lam8s(3)%atype  = 'deg'

! 4. Schicht

lam8s(4)%lamid  = 1
lam8s(4)%plyid  = 4
lam8s(4)%mat8id = 2
lam8s(4)%th	= 0.25D0
lam8s(4)%angle  = 90.D0
lam8s(4)%atype  = 'deg'

! 5. Schicht

lam8s(5)%lamid  = 1
lam8s(5)%plyid  = 5
lam8s(5)%mat8id = 2
lam8s(5)%th	= 0.25D0
lam8s(5)%angle  = 90.D0
lam8s(5)%atype  = 'deg'

! 6. Schicht

lam8s(6)%lamid  = 1
lam8s(6)%plyid  = 6
lam8s(6)%mat8id = 2
lam8s(6)%th	= 0.25D0
lam8s(6)%angle  = 0.D0
lam8s(6)%atype  = 'deg'

! 7. Schicht

lam8s(7)%lamid  = 1
lam8s(7)%plyid  = 7
lam8s(7)%mat8id = 2
lam8s(7)%th	= 0.25D0
lam8s(7)%angle  = 90.D0
lam8s(7)%atype  = 'deg'

! 8. Schicht

lam8s(8)%lamid  = 1
lam8s(8)%plyid  = 8
lam8s(8)%mat8id = 2
lam8s(8)%th	= 0.25D0
lam8s(8)%angle  = 0.D0
lam8s(8)%atype  = 'deg'

! mat8s

allocate(mat8s(1))

mat8s(1)%mid   = 2
mat8s(1)%ym11  = 141000
mat8s(1)%ym22  = 9400
mat8s(1)%nu12  = 0.3
mat8s(1)%sm12  = 4600
mat8s(1)%sm13  = mat8s(1)%sm12
mat8s(1)%sm23  = 3200
mat8s(1)%rho   = 0.D0
mat8s(1)%ath11 = 0.D0
mat8s(1)%ath22 = 0.D0
mat8s(1)%tref  = 0.D0
mat8s(1)%ge    = 0.D0

 ! spcadds

 allocate(spcadds(size(spc1s,1)))

do ii=1,size(spc1s,1)
 
   spcadds(ii)%scid = 1
   spcadds(ii)%sid  = ii

end do
!
!! spc1s
!
! countrb=1

!allocate(spc1s(size(spc1,1)))

! subcases

allocate(subcases(1))

!	 subcases(1) = (/ 1, 1, 1 /)
	 
subcases(1)%scid     = 1
subcases(1)%spcaddid = 1
subcases(1)%loadid   = 1

! solution 1-static, 2-static+buckling

!write(*,*) 'spc1s: ', spc1s
!write(*,*) 'forces: ', forces
!stop

! set estimate for nonzero entries in 1 row of partial stiffness matrices

nzKaaS = 125
nzKabS = 50
nzKbbS = 30

nzKgaaS = 75
!nzKgabS = 5
!nzKgbbS = 10

sol=2
numEigVal = 1

!
! =================================================================================================
!
! Error handling
!

9999 continue

if (err_code /= 0) then
   
   write(*,*)                      'An error occured in subroutine'
   write(*,*)  			   'testcase12'
   write(*,'(A,I2)',advance='YES') ' Errorcode: ', err_code
   write(*,*)                      'exit program '
   stop
   
end if

end subroutine testcase12
