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
subroutine testcase10 (testcase)
			
! =================================================================================================
!
!	Header:		subroutine defining testcase9
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
!	Author:		Martin RÃ¤del			08.12.2009
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
integer, parameter		:: anzelemlaenge=100
integer, parameter		:: anzelembreite=100
double precision, parameter	:: elemlaenge = 1.D0
double precision, parameter	:: elembreite = 1.D0
double precision, parameter	:: fakelemdicke=0.15d0
double precision		:: force=1000.D0

double precision		:: elemdicke
integer, parameter		:: anzelem=anzelemlaenge*anzelembreite
integer, parameter		:: anzknotenbreite=2*anzelembreite+1
integer, parameter		:: anzknotenlaenge=2*anzelemlaenge+1
integer, parameter		:: anzknoten=anzknotenlaenge*anzknotenbreite-anzelem

integer				:: anzspc1

! Krag

!integer, parameter		:: rbx0=0, rbxlaenge=0  						 ! RB an den Seiten (ohne Eckknoten)
!integer, parameter		:: rby0=12345, rbybreite=0
!integer, parameter		:: rbx0y0=12345, rbxlaengey0=12345, rbx0ybreite=0, rbxlaengeybreite=0	 ! Eckknoten
!integer				:: countrb, temp, facI, count_dof

! SS-SS

integer, parameter		:: rbx0=3, rbxlaenge=3  						 ! RB an den Seiten (ohne Eckknoten)
integer, parameter		:: rby0=23, rbybreite=3
integer, parameter		:: rbx0y0=123, rbxlaengey0=23, rbx0ybreite=3, rbxlaengeybreite=3	 ! Eckknoten
integer 			:: countrb, temp, facI, count_dof

! SC-SC

!integer, parameter		 :: rbx0=3, rbxlaenge=35						  ! RB an den Seiten (ohne Eckknoten)
!integer, parameter		 :: rby0=23, rbybreite=34
!integer, parameter		 :: rbx0y0=123, rbxlaengey0=235, rbx0ybreite=34, rbxlaengeybreite=345	  ! Eckknoten
!integer			 :: countrb, temp, facI, count_dof

! SC-SF

!integer, parameter		 :: rbx0=3, rbxlaenge=35						  ! RB an den Seiten (ohne Eckknoten)
!integer, parameter		 :: rby0=23, rbybreite=0
!integer, parameter		 :: rbx0y0=123, rbxlaengey0=235, rbx0ybreite=3, rbxlaengeybreite=35	  ! Eckknoten
!integer			 :: countrb, temp, facI, count_dof

integer				:: num
integer				:: ii,jj

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
if (anzelemlaenge*elemlaenge < anzelembreite*elembreite) then
   elemdicke=fakelemdicke*(anzelemlaenge*elemlaenge)
else
   elemdicke=fakelemdicke*(anzelembreite*elembreite)
end if

write(*,*) 'elemdicke ', elemdicke

!anzf = 6*((anzknotenlaenge-2)+(anzknotenbreite-2))+12
!
! =================================================================================================
!
! Calculation
!
! number of bc's

anzspc1=0

if (rbx0y0 > 0) then
   if (rbx0y0 > 10E+03) then
      count_dof = 5
   else if (rbx0y0 > 10E+02) then
      count_dof = 4
   else if (rbx0y0 > 10E+01) then
      count_dof = 3
   else if (rbx0y0 > 10E+00) then
      count_dof = 2
   else if (rbx0y0 > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+1
end if

if (rbx0 > 0) then
   if (rbx0 > 10E+03) then
      count_dof = 5
   else if (rbx0 > 10E+02) then
      count_dof = 4
   else if (rbx0 > 10E+01) then
      count_dof = 3
   else if (rbx0 > 10E+00) then
      count_dof = 2
   else if (rbx0 > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+(anzknotenbreite-2)
end if

if (rbx0ybreite > 0) then
   if (rbx0ybreite > 10E+03) then
      count_dof = 5
   else if (rbx0ybreite > 10E+02) then
      count_dof = 4
   else if (rbx0ybreite > 10E+01) then
      count_dof = 3
   else if (rbx0ybreite > 10E+00) then
      count_dof = 2
   else if (rbx0ybreite > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+1
end if

if (rbybreite > 0) then
   if (rbybreite > 10E+03) then
      count_dof = 5
   else if (rbybreite > 10E+02) then
      count_dof = 4
   else if (rbybreite > 10E+01) then
      count_dof = 3
   else if (rbybreite > 10E+00) then
      count_dof = 2
   else if (rbybreite > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+1
end if

if (rbxlaengeybreite > 0) then
   if (rbxlaengeybreite > 10E+03) then
      count_dof = 5
   else if (rbxlaengeybreite > 10E+02) then
      count_dof = 4
   else if (rbxlaengeybreite > 10E+01) then
      count_dof = 3
   else if (rbxlaengeybreite > 10E+00) then
      count_dof = 2
   else if (rbxlaengeybreite > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+1
end if

if (rbxlaenge > 0) then
   if (rbxlaenge > 10E+03) then
      count_dof = 5
   else if (rbxlaenge > 10E+02) then
      count_dof = 4
   else if (rbxlaenge > 10E+01) then
      count_dof = 3
   else if (rbxlaenge > 10E+00) then
      count_dof = 2
   else if (rbxlaenge > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+(anzknotenbreite-2)
end if

if (rbxlaengey0 > 0) then
   if (rbxlaengey0 > 10E+03) then
      count_dof = 5
   else if (rbxlaengey0 > 10E+02) then
      count_dof = 4
   else if (rbxlaengey0 > 10E+01) then
      count_dof = 3
   else if (rbxlaengey0 > 10E+00) then
      count_dof = 2
   else if (rbxlaengey0 > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+1
end if

if (rby0 > 0) then
   if (rby0 > 10E+03) then
      count_dof = 5
   else if (rby0 > 10E+02) then
      count_dof = 4
   else if (rby0 > 10E+01) then
      count_dof = 3
   else if (rby0 > 10E+00) then
      count_dof = 2
   else if (rby0 > 0) then
      count_dof = 1
   end if
   anzspc1=anzspc1+1
end if

write(*,*) 'anzspc1: ', anzspc1
!stop
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

! nodes
	 
allocate(nodes(anzknoten))

num = 1

do ii=1,anzknotenbreite
   do jj=1,anzknotenlaenge	! zuerst über länge
      
      if (mod(ii,2) > 0) then 		! ii ungerade
	 
	 nodes(num)%nid       = num
	 nodes(num)%cid       = 0
	 nodes(num)%coords(1) = (jj-1)*elemlaenge/2.D0
	 nodes(num)%coords(2) = (ii-1)*elembreite/2.D0
	 nodes(num)%coords(3) = 0.D0
	 
	 num = num+1
	 
      else
         
	 if (mod(jj,2) > 0) then
	    
	    nodes(num)%nid       = num
	    nodes(num)%cid       = 0
	    nodes(num)%coords(1) = (jj-1)*elemlaenge/2.D0
	    nodes(num)%coords(2) = (ii-1)*elembreite/2.D0
	    nodes(num)%coords(3) = 0.D0
	    
	    num = num+1
	    
	 end if
	 
      end if
      
   end do
end do

! quad8 elements

allocate(quad8s(anzelem))

num=1

do ii=1,anzelembreite
   do jj=1,anzelemlaenge
      
      quad8s(num)%eid     = num
      quad8s(num)%pid     = 1
      quad8s(num)%nids(1) = (jj-1)*3+2-jj+(ii-1)*(anzknotenlaenge+anzelemlaenge+1)
      quad8s(num)%nids(2) = quad8s(num)%nids(1)+2
      quad8s(num)%nids(3) = (jj-1)*3+2-jj+(ii)*(anzknotenlaenge+anzelemlaenge+1)+2
      quad8s(num)%nids(4) = quad8s(num)%nids(3)-2
      quad8s(num)%nids(5) = quad8s(num)%nids(1)+1
      quad8s(num)%nids(6) = quad8s(num)%nids(2)+anzknotenlaenge-jj
      quad8s(num)%nids(7) = quad8s(num)%nids(3)-1
      quad8s(num)%nids(8) = quad8s(num)%nids(1)+anzknotenlaenge-jj+1
      quad8s(num)%theta   = 0.D0
      quad8s(num)%atype   = 'deg'
      quad8s(num)%offset  = 0.D0
      
      num = num+1
      
   end do
end do

!write(*,*) 'quad8s: ', quad8s
!stop
	 
! loads
	 
allocate(loads(anzknotenlaenge))

do ii=1,anzknotenlaenge
	 
loads(ii)%lcid  = 1
loads(ii)%sfac  = 1.D0
loads(ii)%sfaci = 1.D0
loads(ii)%lidi  = ii

end do
	 
! forces
	 
allocate(forces(anzknotenlaenge))

do ii=1,anzknotenlaenge
  
   forces(ii)%lid = ii
   forces(ii)%nid = anzknoten-(ii-1)
   forces(ii)%cid = 0
   if (ii==1 .OR. ii==anzknotenlaenge) then
      forces(ii)%fac = 1.D0/6.D0*force
   else
      
      if (mod(ii,2) > 0) then
	 forces(ii)%fac = 1.D0/3.D0*force
      else
         forces(ii)%fac = 2.D0/3.D0*force
      end if
      
   end if
   forces(ii)%ni  = (/ 0.D0, 1.D0, 0.D0 /)
  
end do
	 
! pshells

allocate(pshells(1))
        
pshells(1)%pid  = 1
pshells(1)%mid1 = 1
pshells(1)%mt	= elemdicke
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
	 
mat1s(1)%mid  = 1
mat1s(1)%ym   = 210000.D0
mat1s(1)%sm   = 78947.37D0
mat1s(1)%nu   = 0.33D0
mat1s(1)%rho  = 0.D0
mat1s(1)%ath  = 0.D0
mat1s(1)%tref = 0.D0
mat1s(1)%ge   = 0.D0

! pcomps

allocate(pcomps(1))

pcomps(1)%pid    = 2
pcomps(1)%lamid  = 1
pcomps(1)%offset = 0.D0
pcomps(1)%lay    = 3
pcomps(1)%nsm    = 0
pcomps(1)%sb     = 0.D0

! lam8s

allocate(lam8s(3))

! 1. Schicht

lam8s(1)%lamid  = 1
lam8s(1)%plyid  = 1
lam8s(1)%mat8id = 2
lam8s(1)%th     = 0.25D0
lam8s(1)%angle  = 0.D0
lam8s(1)%atype  = 'deg'

! 2. Schicht

lam8s(2)%lamid  = 1
lam8s(2)%plyid  = 2
lam8s(2)%mat8id = 2
lam8s(2)%th     = 0.25D0
lam8s(2)%angle  = 90.D0
lam8s(2)%atype  = 'deg'

! 3. Schicht

lam8s(3)%lamid  = 1
lam8s(3)%plyid  = 3
lam8s(3)%mat8id = 2
lam8s(3)%th     = 0.25D0
lam8s(3)%angle  = 45.D0
lam8s(3)%atype  = 'deg'

! 4. Schicht

!lam8s(4)%lamid  = 1
!lam8s(4)%plyid  = 4
!lam8s(4)%mat8id = 2
!lam8s(4)%th     = 0.25D0
!lam8s(4)%angle  = 45.D0
!lam8s(4)%atype  = 'deg'

! 5. Schicht

!lam8s(5)%lamid  = 1
!lam8s(5)%plyid  = 5
!lam8s(5)%mat8id = 2
!lam8s(5)%th	 = 0.25D0
!lam8s(5)%angle  = 90.D0
!lam8s(5)%atype  = 'deg'

! 6. Schicht

!lam8s(6)%lamid  = 1
!lam8s(6)%plyid  = 6
!lam8s(6)%mat8id = 2
!lam8s(6)%th	 = 0.25D0
!lam8s(6)%angle  = 0.D0
!lam8s(6)%atype  = 'deg'

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

allocate(spcadds(anzspc1))

do ii=1,anzspc1
        
   spcadds(ii)%scid = 1
   spcadds(ii)%sid  = ii

end do

! spc1s

countrb=1

allocate(spc1s(anzspc1))

! RB bei x=0, y=0

if (rbx0y0 > 0) then

   spc1s(1)%sid = 1
   spc1s(1)%dof = rbx0y0
   spc1s(1)%n1  = 1
   spc1s(1)%thru= .false.
   spc1s(1)%nn  = spc1s(1)%n1

   countrb=countrb+1

end if

! RB bei x=0 (außer Eckknoten)

if (rbx0 > 0) then

   temp=countrb
   
   jj = 2
   num = 1+anzknotenlaenge

   do ii=countrb,(countrb+anzknotenbreite-3)

      spc1s(ii)%sid = 1
      spc1s(ii)%dof = rbx0
      spc1s(ii)%n1  = num
      spc1s(ii)%thru= .false.
      spc1s(ii)%nn  = spc1s(ii)%n1
   
      if (mod(jj,2) > 0) then
      
	 num = num+anzknotenlaenge
	 
      else
	 
	 num = num+anzelemlaenge+1
	 
      end if

      temp=temp+1
      jj = jj+1
      
   end do

   countrb=temp
 
end if

! RB bei x=0, y=breite

if (rbx0ybreite > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rbx0ybreite
   spc1s(countrb)%n1  = anzknoten-anzknotenlaenge+1
   spc1s(countrb)%thru= .false.
   spc1s(countrb)%nn  = spc1s(countrb)%n1

    countrb=countrb+1
    
end if
    
! RB bei y=breite (außer Eckknoten)
   
if (rbybreite > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rbybreite
   spc1s(countrb)%n1  = anzknoten-anzknotenlaenge+2
   spc1s(countrb)%thru= .true.
   spc1s(countrb)%nn  = anzknoten-1

    countrb=countrb+1

end if
    
! RB bei x=laenge, y=breite

if (rbxlaengeybreite > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rbxlaengeybreite
   spc1s(countrb)%n1  = anzknoten
   spc1s(countrb)%thru= .false.
   spc1s(countrb)%nn  = spc1s(countrb)%n1

    countrb=countrb+1

end if
    
! RB bei x=laenge (außer Eckknoten)

if (rbxlaenge > 0) then

   temp=countrb
   jj=2
   num = anzknoten-anzknotenlaenge
   
   do ii=countrb,(countrb+anzknotenbreite-3)
   
      spc1s(ii)%sid = 1
      spc1s(ii)%dof = rbxlaenge
      spc1s(ii)%n1  = num
      spc1s(ii)%thru= .false.
      spc1s(ii)%nn  = spc1s(ii)%n1
      
      if (mod(jj,2) > 0) then
      
	 num = num-anzknotenlaenge
	 
      else
	 
	 num = num-(anzelemlaenge+1)
	 
      end if
   
      temp=temp+1
      jj = jj+1
      
   end do

   countrb=temp

end if

! RB bei x=laenge, y=0

if (rbxlaengey0 > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rbxlaengey0
   spc1s(countrb)%n1  = anzknotenlaenge
   spc1s(countrb)%thru= .false.
   spc1s(countrb)%nn  = spc1s(countrb)%n1

    countrb=countrb+1 

end if

! RB bei y=0 (außer Eckknoten)

if (rby0 > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rby0
   spc1s(countrb)%n1  = 2
   spc1s(countrb)%thru= .true.
   spc1s(countrb)%nn  = anzknotenlaenge-1

end if

!! spcadds

!allocate(spcadds(2))

!do ii=1,2
	 
!   spcadds(ii)%scid = 1
!   spcadds(ii)%sid  = ii

!end do

!allocate(spc1s(1))

!spc1s(1)%sid = 1
!spc1s(1)%dof = 12345
!spc1s(1)%n1  = 1
!spc1s(1)%thru= .true.
!spc1s(1)%nn  = 3

!spc1s(2)%sid = 1
!spc1s(2)%dof = 12345
!spc1s(2)%n1  = 3
!spc1s(2)%thru= .false.
!spc1s(2)%nn  = 3

if (countrb == anzspc1) then
      
   write(*,*) 'Modellaufbau fertig'
   
else
   
   write(*,*) 'Fehler beim Modellaufbau'
   write(*,*) 'Aufbau der spc1-Karte fehlgeschlagen'
   write(*,*) 'spc1s', spc1s
   !stop
   
end if

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

nzKaaS = 150
nzKabS = 75
nzKbbS = 50

nzKgaaS = 75
!nzKgabS = 5
!nzKgbbS = 10

sol=2
numEigVal = 3

end subroutine testcase10
