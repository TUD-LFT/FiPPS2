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
subroutine testcase8 (testcase)
			
! =================================================================================================
!
!	Header:		subroutine defining testcase4
!
!	Content:	
!			 tria3-element in global coordinates
!			 force in -z-direction on node 
!			 x,y,z,rotx,roty-bc on nodes 1,2
!
!			slender cantilever
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

! inner


integer				:: ii,jj,nn
integer,parameter		:: anzelemprohalbbreite=50
integer,parameter		:: anzelemprohalblaenge=50

double precision,parameter	:: fakelemdicke=0.2d0
double precision		:: elemdicke

integer,parameter		:: anzelem=8*anzelemprohalblaenge*anzelemprohalbbreite

integer,parameter		:: anzknotenbreite=2*anzelemprohalbbreite+1
integer,parameter		:: anzknotenlaenge=2*anzelemprohalblaenge+1
integer,parameter		:: anzknoten=anzknotenlaenge*anzknotenbreite
!integer,parameter		:: anzspc1=6+2*(anzknotenbreite-2)
integer				:: anzspc1
integer,parameter		:: rbx0=3, rbxlaenge=3						! RB an den Seiten (ohne Eckknoten)
integer,parameter		:: rby0=23, rbybreite=3
integer,parameter		:: rbx0y0=123, rbxlaengey0=23, rbx0ybreite=3, rbxlaengeybreite=3	! Eckknoten
integer				:: countrb, temp, facI, count_dof
!
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
if (anzelemprohalbbreite < anzelemprohalblaenge) then
   elemdicke=fakelemdicke*(anzelemprohalbbreite*2)
else
   elemdicke=fakelemdicke*(anzelemprohalblaenge*2)
end if

write(*,*) 'elemdicke ', elemdicke
!stop
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
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
   anzspc1=anzspc1+count_dof*(anzknotenbreite-2)
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
   anzspc1=anzspc1+count_dof
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
   anzspc1=anzspc1+count_dof*(anzknotenbreite-2)
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
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
!   write(*,*) 'anzspc1: ', anzspc1
!   write(*,*) 'count_dof: ', count_dof
!   stop
end if

write(*,*) 'anzspc1: ', anzspc1
!stop


!
! =================================================================================================

! 
         
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
	 
allocate(nodes(anzknoten))

do ii=1,anzknotenbreite
   do jj=1,anzknotenlaenge
      
      nn = (ii-1)*anzknotenlaenge+jj
   
      nodes(nn)%nid       = nn
      nodes(nn)%cid       = 0
      nodes(nn)%coords(1) = DBLE(jj-1)
      nodes(nn)%coords(2) = DBLE(ii-1)
      nodes(nn)%coords(3) = 0
      
   end do
   
end do

! tria3 elements

allocate(tria3s(anzelem))

do ii=1,anzknotenlaenge-1
   do jj=1,anzknotenbreite-1
      
      nn = 2*((ii-1)*(anzknotenbreite-1)+jj)-1
      
      if ((ii >  anzelemprohalblaenge .AND. jj <= anzelemprohalbbreite) .OR. &
      	& (ii <= anzelemprohalblaenge .AND. jj >  anzelemprohalbbreite)) then
	 
	    tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) = (jj-1)*anzknotenlaenge+ii
            tria3s(nn)%nids(2) = (jj-1)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(3) =   (jj)*anzknotenlaenge+ii+1
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%atype   = 'deg'
	    tria3s(nn)%offset  = 0.D0
            
            nn=nn+1
            
            tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) =   (jj)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(2) =   (jj)*anzknotenlaenge+ii
            tria3s(nn)%nids(3) = (jj-1)*anzknotenlaenge+ii
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%atype   = 'deg'
	    tria3s(nn)%offset  = 0.D0
	
      else
   
            tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) = (jj-1)*anzknotenlaenge+ii
            tria3s(nn)%nids(2) = (jj-1)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(3) =   (jj)*anzknotenlaenge+ii
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%atype   = 'deg'
	    tria3s(nn)%offset  = 0.D0
            
            nn=nn+1
            
            tria3s(nn)%eid     = nn
            tria3s(nn)%pid     = 1
            tria3s(nn)%nids(1) =   (jj)*anzknotenlaenge+ii+1
            tria3s(nn)%nids(2) =   (jj)*anzknotenlaenge+ii
            tria3s(nn)%nids(3) = (jj-1)*anzknotenlaenge+ii+1
	    tria3s(nn)%theta   = 0.D0
	    tria3s(nn)%atype   = 'deg'
	    tria3s(nn)%offset  = 0.D0
	    
      end if
            
    end do
   
end do
	 
! loads
	 
!allocate(loads(anzknotenlaenge))

!do ii=1,anzknotenlaenge
	 
!loads(ii)%lcid  = 1
!loads(ii)%sfac  = 1.D0
!loads(ii)%sfaci = 1.D0
!loads(ii)%lidi  = ii

!end do

! loads
	 
allocate(loads(1))
	 
loads(1)%lcid  = 1
loads(1)%sfac  = 1.D0
loads(1)%sfaci = 1.D0
loads(1)%lidi  = 1
	 
! forces
	 
!allocate(forces(anzknotenlaenge))

!do ii=1,anzknotenlaenge
   
!   forces(ii)%lid = 1
!   forces(ii)%nid = anzknoten-(ii-1)
!   forces(ii)%cid = 0
!   if (ii==1 .OR. ii==anzknotenlaenge) then
!      forces(ii)%fac = 5.D0
!   else
!      forces(ii)%fac = 10.D0
!   end if
!   forces(ii)%ni  = (/ 0.D0, -2.D0, 0.D0 /)
   
!end do

! forces
	 
allocate(forces(anzknotenlaenge))
   
do ii=1,anzknotenlaenge
  
   forces(ii)%lid = 1
   forces(ii)%nid = anzknoten/2+anzknotenlaenge/2-(ii-1)
   forces(ii)%cid = 0
   if (ii==1 .OR. ii==anzknotenlaenge) then
      forces(ii)%fac = 5.D0
   else
      forces(ii)%fac = 10.D0
   end if
   forces(ii)%ni  = (/ 0.D0, 0.D0, -1000.D0 /)
  
end do
	 
! pshells

allocate(pshells(1))
	 
!	 pshells(1) = (/ 1, 1, 1.D0, 1, 1.D0, 1, 0.833333333D0, 0.D0, -0.5D0, 0.5D0, 1 /)
	 
pshells(1)%pid  = 1
pshells(1)%mid1 = 1
pshells(1)%mt   = elemdicke
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

!allocate(spcadds(anzknotenbreite))
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

   do ii=countrb,(countrb+anzknotenbreite-3)

      spc1s(ii)%sid = 1
      spc1s(ii)%dof = rbx0
      spc1s(ii)%n1  = (ii-1)*anzknotenlaenge+1
      spc1s(ii)%thru= .false.
      spc1s(ii)%nn  = spc1s(ii)%n1
   !   spc1s(ii)%inc = 1

      temp=temp+1
      
   end do

   countrb=temp
 
end if

! RB bei x=0, y=breite

if (rbx0ybreite > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rbx0ybreite
!   spc1s(countrb)%n1  = (anzknotenbreite-1)*anzknotenlaenge+1
   spc1s(countrb)%n1  = anzknoten-anzknotenlaenge+1
   spc1s(countrb)%thru= .false.
   spc1s(countrb)%nn  = spc1s(countrb)%n1

    countrb=countrb+1
    
end if
    
! RB bei y=breite (außer Eckknoten)
   
if (rbybreite > 0) then

   spc1s(countrb)%sid = 1
   spc1s(countrb)%dof = rbybreite
   spc1s(countrb)%n1  = (anzknotenbreite-1)*anzknotenlaenge+2
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
   facI=1
   
   do ii=countrb,(countrb+anzknotenbreite-3)
   
      spc1s(ii)%sid = 1
      spc1s(ii)%dof = rbxlaenge
      spc1s(ii)%n1  = anzknoten-facI*anzknotenlaenge
      spc1s(ii)%thru= .false.
      spc1s(ii)%nn  = spc1s(ii)%n1
   !   spc1s(ii)%inc = 1
   
      temp=temp+1
      facI=facI+1
      
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

if (countrb == anzspc1) then
      
   write(*,*) 'Modellaufbau fertig'
!   write(*,*) 'spc1s', spc1s
!   write(*,*) 'anzspc1', anzspc1
!   write(*,*) 'countrb', countrb
!   stop
   
else
   
   write(*,*) 'Fehler beim Modellaufbau'
   write(*,*) 'Aufbau der spc1-Karte fehlgeschlagen'
   write(*,*) 'spc1s', spc1s
   stop
   
end if

! subcases

allocate(subcases(1))

!	 subcases(1) = (/ 1, 1, 1 /)
	 
subcases(1)%scid     = 1
subcases(1)%spcaddid = 1
subcases(1)%loadid   = 1

! solution 1-static, 2-static+buckling

nzKaaS = 40
nzKabS = 10
nzKbbS = 15

nzKgaaS = 20

sol=1
numEigVal = 10

end subroutine testcase8
